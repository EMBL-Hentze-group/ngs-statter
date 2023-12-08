FROM python:3.11-slim-bookworm
COPY requirements.txt .
ARG USERNAME=sahadeva
ARG USER_UID=1001
ARG USER_GID=$USER_UID
ENV CARGO_HOME=/usr/local/cargo
ENV RUSTUP_HOME=/usr/local/rustup
RUN apt-get update \
    && apt-get install cmake curl libfindbin-libs-perl -y \
    && curl https://sh.rustup.rs -sSf |  sh -s -- --default-toolchain stable -y \
    && pip install --no-cache-dir -r requirements.txt \
    && groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && chown -R sahadeva:sahadeva  /usr/local/cargo /usr/local/cargo
ENV PATH=/usr/local/cargo/bin:$PATH
USER $USERNAME