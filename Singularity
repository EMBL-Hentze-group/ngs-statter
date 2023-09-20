Bootstrap: docker
From: python:3.11-slim-bookworm
IncludeCmd: yes

%environment
    LC_ALL="en_US.UTF-8"
    export LC_ALL

%labels
    Author Sudeep Sahadevan
    build_date 2023 September 11

%runscript
    exec "$@"

%files
    statter/* statter/
%post
    cd statter
    pip install .
    rm -rf build dist mad_statter.egg-info
    pip cache purge