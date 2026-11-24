#!/bin/sh

# This script installs files [$3 ...] into the directory [$1] with
# the install mode [$2].  Each file keeps its path relative to the
# current directory, so src/decays/decay.cpp is installed as
# [$1]/src/decays/decay.cpp.  Empty arguments are ignored.

if test $# -lt 2; then
    echo "Usage: ./`basename $0` destination mode [files]"
    exit 1
fi

dest="$1"
shift
mode="$1"
shift

for file in "$@"; do
    case "${file}" in
        "") continue ;;
        /*) echo "Error: absolute path not supported: ${file}"
            exit 1 ;;
    esac

    if test ! -e "${file}"; then
        echo "Error: file does not exist: ${file}"
        exit 1
    fi

    dir="${dest}/`dirname ${file}`"

    install -d "${dir}" || exit 1
    install -m "${mode}" "${file}" "${dir}" || exit 1
done
