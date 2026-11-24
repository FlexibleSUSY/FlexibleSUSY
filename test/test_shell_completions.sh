#!/bin/bash

# Checks that the bash and zsh completions know the options of configure,
# createmodel, createaddon and the spectrum generators.

BASEDIR=$(dirname "$0")
ROOT="${BASEDIR}/.."
BASH_COMPLETION="${ROOT}/utils/install-bash_completions.bash"
ZSH_COMPLETION="${ROOT}/utils/install-zsh_completions.zsh"

errors=0

error() {
    echo "Error: $*"
    errors=$((errors + 1))
}

# checks that a file contains an option, delimited by something that
# can't be part of an option name
check_contains() {
    local file="$1" option="$2"
    grep -Eq -- "(^|[^-a-z0-9])${option}([^-a-z0-9]|$)" "${file}" ||
        error "${option} is missing in $(basename "${file}")"
}

# checks that all --with-* options of a completion are known
check_known() {
    local file="$1" known="$2" option
    for option in $(grep -oE -- '--with-[a-z0-9-]+' "${file}" | sort -u); do
        grep -qx -- "${option}" <<<"${known}" ||
            error "${option} in $(basename "${file}") is not a configure option"
    done
}

configure_help=$("${ROOT}/configure" --help)

# --with-* options of configure
with_options=$(grep -oE '^ +--with-[a-z0-9-]+=' <<<"${configure_help}" |
    sed -E 's/^ +//; s/=$//' | sort -u)

# enable/disable features of configure
features=$(sed -n '/^enable\/disable options/,/^Package directories/p' <<<"${configure_help}" |
    grep -oE '^  [a-z][a-z0-9-]+' | sed 's/^ *//')

[ -n "${with_options}" ] || error "could not read the --with-* options of configure"
[ -n "${features}" ] || error "could not read the features of configure"

for option in ${with_options}; do
    check_contains "${BASH_COMPLETION}" "${option}"
    check_contains "${ZSH_COMPLETION}" "${option}"
done

for feature in ${features}; do
    check_contains "${BASH_COMPLETION}" "--enable-${feature}"
    check_contains "${BASH_COMPLETION}" "--disable-${feature}"
    # the zsh script creates --enable-* and --disable-* from this list
    grep -q "^ *'${feature}\[" "${ZSH_COMPLETION}" ||
        error "${feature} is missing in $(basename "${ZSH_COMPLETION}")"
done

check_known "${BASH_COMPLETION}" "${with_options}"
check_known "${ZSH_COMPLETION}" "${with_options}"

# options of createmodel and createaddon
for script in createmodel createaddon; do
    for option in $("${ROOT}/${script}" --help |
        grep -oE '^  --[a-z-]+' | sed 's/^ *//' | sort -u); do
        check_contains "${BASH_COMPLETION}" "${option}"
        check_contains "${ZSH_COMPLETION}" "${option}"
    done
done

# options of the spectrum generators
for option in $(grep -oE '"  --[a-z-]+' "${ROOT}/src/command_line_options.cpp" |
    sed 's/^"  //' | sort -u); do
    check_contains "${BASH_COMPLETION}" "${option}"
    check_contains "${ZSH_COMPLETION}" "${option}"
done

bash -n "${BASH_COMPLETION}" || error "syntax error in the bash completion"

if command -v zsh >/dev/null 2>&1; then
    zsh -n "${ZSH_COMPLETION}" || error "syntax error in the zsh completion"
fi

if [ "${errors}" -ne 0 ]; then
    echo "Test FAILED: ${errors} error(s)"
    exit 1
fi

echo "All tests passed"
