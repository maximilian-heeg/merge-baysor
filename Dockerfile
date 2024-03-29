FROM rust AS build
WORKDIR /app
COPY . .
RUN cargo build --release 
RUN strip target/release/merge-baysor
ENTRYPOINT ["/app/target/release/merge-baysor"]
CMD ["--help"]


FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY --from=build /app/target/release/merge-baysor /app/
ENV PATH="${PATH}:/app/"
ENTRYPOINT ["/app/merge-baysor"]
CMD ["--help"]