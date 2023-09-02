FROM rust AS build
WORKDIR /app
COPY . .
RUN cargo build --release 
RUN strip target/release/merge-baysor
ENTRYPOINT ["/app/target/release/merge-baysor"]
CMD ["--help"]


FROM debian:bookworm-slim
WORKDIR /app
COPY --from=build /app/target/release/merge-baysor /app/
ENTRYPOINT ["/app/merge-baysor"]
CMD ["--help"]