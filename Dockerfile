FROM python:3.12-slim-trixie AS build
COPY . /mad_statter
ENV CARGO_HOME=/usr/local/cargo
ENV RUSTUP_HOME=/usr/local/rustup
ENV PATH=/usr/local/cargo/bin:$PATH
RUN apt-get update \
    && apt-get install cmake curl libfindbin-libs-perl -y \
    && curl https://sh.rustup.rs -sSf |  sh -s -- --default-toolchain stable -y \
    && pip install --no-cache-dir poetry \
    && cd mad_statter \
    && poetry config virtualenvs.create false \
    && poetry install  \
    && poetry build

FROM python:3.12-slim-trixie
COPY --from=build /mad_statter/dist/* /tmp/
RUN cd /tmp/ \
    && ls *.whl| xargs -I whl pip install --no-cache-dir whl \
    && rm -rf *.whl \
    && pip cache purge
