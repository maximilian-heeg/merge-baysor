use clap::Parser;
use polars::prelude::*;
use std::error::Error;
use std::fmt;
use std::fs::File;

/// Merge segmentations results from different runs, e.g. FOV
/// into a single segmetation file. This is done by combining cells if the intersection over union
/// is bigger than a set threshold. 
/// This approach is similar to the stichting that is done in cellpose.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None, verbatim_doc_comment)]
pub struct Args {
    /// The Baysor output files.
    /// These files must include the following columns:
    /// - transcript_id
    /// - cell
    #[arg(verbatim_doc_comment)]
    files: Vec<String>,

    /// Threshold for stitching. If the IOU for two cells is greater that the threshold, they will be merged.
    #[arg(long, default_value_t = 0.2)]
    threshold: f32,

    /// Additional columns that will be included in the final output CSV file.
    #[arg(long, default_values = ["x", "y", "z", "qv", "overlaps_nucleus", "gene"])]
    additional_columns: Vec<String>,

    /// Output file
    #[arg(long, default_value = "out.csv")]
    outfile: String,
}


#[derive(Debug)]
struct MyError(String);

impl fmt::Display for MyError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "There is an error: {}", self.0)
    }
}
impl Error for MyError {}


/// The main logic of that program lies here.
/// It reads in the files, reduces the list by merging the layers, and the combines it with the additional columns to create a output file.
pub fn run(args: Args) -> Result<(), Box<dyn Error>> {
    if args.files.len() < 2 {
        return Err(Box::new(MyError(
            "Not enought files to merge. Please provide at least two files.".into(),
        )));
    }

    // Read in all the layers
    // TODO: Improvement. Make sure that cell_id are unique. They should be Baysor created cell_ids based in the uuid of the process.
    println!("Read files");
    let layers: Vec<LazyFrame> = args
        .files
        .iter()
        .map(|v| {
            LazyCsvReader::new(v)
                .has_header(true)
                .finish()
                .unwrap()
                .select(&[col("transcript_id"), col("cell")])
        })
        .collect();

    // The main part of the script. Reduce the Vector by calling our merge function.
    println!("Merging files");
    let df_result = layers.into_iter().reduce(|a, b| merge(a, b, &args));
    let df_result = match df_result {
        Some(x) => x,
        None => {
            return Err(Box::new(MyError(
                "Error occured while merging the layers".into(),
            )))
        }
    };

    // Create a list of all unique transcripts with the additional columns.
    let df_all_transcripts = unique_transcripts(&args.files, &args.additional_columns)?.left_join(
        df_result,
        col("transcript_id"),
        col("transcript_id"),
    );
    let mut df = df_all_transcripts.collect()?;

    println!("All done....");
    println!("{:?}", df);

    println!("Saving to \"{}\"", args.outfile);
    let mut output_file: File = File::create(args.outfile).unwrap();
    CsvWriter::new(&mut output_file).finish(&mut df).unwrap();

    Ok(())
}

/// Merge two layers of transcripts
/// For all cell pairs create a IOU. If a cell is new, keep it. If a cell has an IOU > threshold, merge it.
fn merge(lhs: LazyFrame, rhs: LazyFrame, args: &Args) -> LazyFrame {
    // Create a full table of all transcripts from the two layers: LHS and RHS
    let df_join = lhs.outer_join(rhs, col("transcript_id"), col("transcript_id"));

    // Create transcript counts per cell for both RHS and LHS cells
    // This information will be used to calculate the IOU
    let df_tally = df_join
        .clone()
        .groupby(["cell", "cell_right"])
        .agg([col("transcript_id").count().alias("transcript_counts")]);
    let df_tally_a = df_tally
        .clone()
        .groupby(["cell"])
        .agg([col("transcript_counts").sum().alias("tally_a")]);
    let df_tally_b = df_tally
        .clone()
        .groupby(["cell_right"])
        .agg([col("transcript_counts").sum().alias("tally_b")]);

    // IOU = (transcripts in A and B) / (transcrips in A + transcripts in B - transcripts in A and B)
    // where A is the cell on the LHS and B is the cell on the RHS.
    let df_iou = df_tally
        // join the total counts per cell for each of the two tables
        .left_join(df_tally_a, col("cell"), col("cell"))
        .left_join(df_tally_b, col("cell_right"), col("cell_right"))
        // calcuate the IOU as in cellpose
        .with_columns([(col("transcript_counts").cast(DataType::Float32)
            / (col("tally_a") + col("tally_b")
                - col("transcript_counts").cast(DataType::Float32)))
        .alias("iou")]);

    // Find cells that are new. They sould have the biggest IOU with the null cell id on the lhs
    let df_new_cells = find_new_cells(df_iou.clone());

    // Find cells that should be merged.
    let df_merge_cells = find_cells_to_merge(df_iou, args.threshold);

    // Create the final results table.
    let df_result = df_join
        .left_join(
            df_new_cells
                .clone()
                .select(&[col("cell_right"), col("rank_cell_right").alias("is_new")]),
            col("cell_right"),
            col("cell_right"),
        )
        .left_join(
            df_merge_cells
                .clone()
                .select(&[col("cell_right"), col("cell").alias("merged_cell_id")]),
            col("cell_right"),
            col("cell_right"),
        )
        .with_columns([when(col("merged_cell_id").is_not_null())
            .then(col("merged_cell_id"))
            .otherwise(
                when(col("is_new").eq(1))
                    .then(col("cell_right"))
                    .otherwise(col("cell")),
            )
            .alias("new_cell")])
        .select(&[col("transcript_id"), col("new_cell").alias("cell")]);
    println!(".. working");
    // get the results for the first merge.. then move on to the next round...
    let df_result = df_result.collect().unwrap();
    df_result.lazy()
}

/// Takes a Lazy DataFrame and finds cells, that are only in the rhs but not in the lhs.
/// For this, cell_right is filtered on non-null cells, that ranked by IOU.
/// If a cell is not present on the lhs, the hightest IOU should be with cell "null" on the lhs
fn find_new_cells(df_iou: LazyFrame) -> LazyFrame {
    df_iou
        .filter(col("cell_right").is_not_null())
        .select(&[col("cell"), col("cell_right"), col("iou")])
        // Create a rank for hightst matching cell for both cells (left and right)
        .with_columns([col("iou")
            .rank(
                RankOptions {
                    method: RankMethod::Max,
                    descending: true,
                },
                Some(0),
            )
            .over(&[col("cell_right")])
            .alias("rank_cell_right")])
        .filter(col("rank_cell_right").eq(1).and(col("cell").is_null()))
}

/// Find cells to merge
/// Empty cells on lhs and rhs are removed
/// Ranks are created for the best match for both cells on the lhs and rhs.
/// if the rank==1 for both sides, the cells will be merged
fn find_cells_to_merge(df_iou: LazyFrame, threshold: f32) -> LazyFrame {
    df_iou
        // Remove null cells, these we not be merged anyway
        .filter(col("cell").is_not_null())
        .filter(col("cell_right").is_not_null())
        // Subset by threshold
        .filter(col("iou").gt_eq(threshold))
        .select(&[col("cell"), col("cell_right"), col("iou")])
        // Create a rank for hightst matching cell for both cells (left and right)
        .with_columns([
            col("iou")
                .rank(
                    RankOptions {
                        method: RankMethod::Max,
                        descending: true,
                    },
                    Some(0),
                )
                .over(&[col("cell")])
                .alias("rank_cell"),
            col("iou")
                .rank(
                    RankOptions {
                        method: RankMethod::Max,
                        descending: true,
                    },
                    Some(0),
                )
                .over(&[col("cell_right")])
                .alias("rank_cell_right"),
        ])
        // Only keep the best matching cells for both
        .filter(col("rank_cell").eq(1))
        .filter(col("rank_cell_right").eq(1))
}

/// Build a list of unique transcripts from all passed subsets
/// Include additional columns
fn unique_transcripts(files: &Vec<String>, columns: &Vec<String>) -> Result<LazyFrame, Box<dyn Error>> {
    let result: Vec<LazyFrame> = files
        .iter()
        .map(|v| {
            LazyCsvReader::new(v)
                .has_header(true)
                .finish()
                .unwrap()
                .select(&[
                    col("transcript_id"),
                    cols(columns),
                ])
        })
        .collect();

    let df_concat = diag_concat_lf(result, true, true)?.unique(None, UniqueKeepStrategy::First);

    Ok(df_concat)
}
