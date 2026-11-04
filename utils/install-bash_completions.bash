# This script installs bash completions for FlexibleSUSY's configure,
# createmodel, createaddon scripts as well as all spectrum generators.
#
# Usage:
#
#   . install-bash_completions.bash

__find_available_models()
{
    find $(dirname $0)/model_files/ -mindepth 2 -type f -name FlexibleSUSY.m.in -exec bash -c 'basename $(dirname $1)' modelname {} \; | sort -u
}

__find_created_models()
{
    find $(dirname $0)/models/ -mindepth 2 -type f -name module.mk -exec bash -c 'basename $(dirname $1)' modelname {} \; | sort -u
}

__find_available_addons()
{
    find $(dirname $0)/addons/ -mindepth 1 -type d -name "*" -exec basename {} \; | sort -u
}

__build_filename_completion_list()
{
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts="$@"

    [ "x$opts" = x ] && return 0

    # handle --xxxxxx=
    if [[ ${prev} == "--"* && ${cur} == "=" ]] ; then
        compopt -o filenames
        COMPREPLY=(*)
        return 0
    fi

    # handle --xxxxx=path
    if [[ ${prev} == '=' ]] ; then
        # unescape space
        cur=${cur//\\ / }
        # expand ~ to $HOME
        [[ ${cur} == "~/"* ]] && cur=${cur/\~/$HOME}
        # show completion if path exists (and escape spaces)
        compopt -o filenames
        local files=("${cur}"*)
        [[ -e ${files[0]} ]] && COMPREPLY=( "${files[@]// /\ }" )
        return 0
    fi

    # handle other options
    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
    if [[ ${#COMPREPLY[@]} == 1 && ${COMPREPLY[0]} != "--"*"=" ]] ; then
        # if there's only one option, without =, then allow a space
        compopt +o nospace
    fi

    return 0
}

__select_from_list()
{
    local pprev="$1"
    shift
    local prev="$1"
    shift
    local cur="$1"
    shift
    local option="$1"
    shift
    local list="$@"

    # handle --XXX=
    if [[ ${prev}${cur} == "${option}" ]] ; then
        echo "$list"
        return 0
    fi

    # handle --XXX=YYY
    if [[ ${pprev}${prev} == "${option}" ]] ; then
        local x
        local selected_list=$(for x in ${list} ; do [[ $x == "${cur}"* ]] && echo $x; done)
        echo "$selected_list"
        return 0
    fi
}

# same as __select_from_list, for an option that takes a comma separated list:
# completes the last entry and doesn't offer entries that are already selected
__select_from_comma_list()
{
    local pprev="$1"
    shift
    local prev="$1"
    shift
    local cur="$1"
    shift
    local option="$1"
    shift
    local list="$@"

    # handle --XXX=
    if [[ ${prev}${cur} == "${option}" ]] ; then
        echo "$list"
        return 0
    fi

    # handle --XXX=YYY,ZZZ
    if [[ ${pprev}${prev} == "${option}" ]] ; then
        local selected=""
        [[ ${cur} == *,* ]] && selected="${cur%,*},"
        local last="${cur##*,}"
        local x
        for x in ${list} ; do
            [[ ",${selected}" == *",${x},"* ]] && continue
            [[ $x == "${last}"* ]] && echo "${selected}${x}"
        done
        return 0
    fi
}

_run_spectrum_generator()
{
    local opts="
--slha-input-file=
--slha-output-file=
--spectrum-output-file=
--database-output-file=
--higgsbounds-dataset=
--higgssignals-dataset=
--lilith-db=
--rgflow-output-file=
--build-info
--model-info
--help
--version
"

    __build_filename_completion_list "${opts}"
}

_configure()
{
    local cur prev pprev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    pprev="${COMP_WORDS[COMP_CWORD-2]}"
    opts="
--enable-colors
--enable-compile
--enable-compiler-warnings
--enable-debug
--enable-feynarts
--enable-fflite
--enable-formcalc
--enable-gm2calc
--enable-higgstools
--enable-himalaya
--enable-librarylink
--enable-lilith
--enable-looptools
--enable-mass-error-check
--enable-meta
--enable-shared-libs
--enable-silent
--enable-sqlite
--enable-static
--enable-threads
--enable-tsil
--enable-verbose

--disable-colors
--disable-compile
--disable-compiler-warnings
--disable-debug
--disable-feynarts
--disable-fflite
--disable-formcalc
--disable-gm2calc
--disable-higgstools
--disable-himalaya
--disable-librarylink
--disable-lilith
--disable-looptools
--disable-mass-error-check
--disable-meta
--disable-shared-libs
--disable-silent
--disable-sqlite
--disable-static
--disable-threads
--disable-tsil
--disable-verbose

--with-addons=
--with-boost-libdir=
--with-boost-incdir=
--with-collier-libdir=
--with-collier-incdir=
--with-cppflags=
--with-cxx=
--with-cxxflags=
--with-eigen-incdir=
--with-fc=
--with-fflags=
--with-flibs=
--with-gm2calc-libdir=
--with-gm2calc-incdir=
--with-gsl-config=
--with-higgstools-libdir=
--with-higgstools-incdir=
--with-himalaya-libdir=
--with-himalaya-incdir=
--with-install-dir=
--with-lilith=
--with-llflags=
--with-lllibs=
--with-loop-libraries=
--with-looptools-libdir=
--with-looptools-incdir=
--with-math-cmd=
--with-models=
--with-optional-modules=
--with-pthread-libdir=
--with-shared-ldflags=
--with-shared-ldlibs=
--with-shared-lib-cmd=
--with-shared-lib-ext=
--with-solvers=
--with-sqlite-libdir=
--with-sqlite-incdir=
--with-static-ldflags=
--with-static-ldlibs=
--with-static-lib-cmd=
--with-static-lib-ext=
--with-tsil-libdir=
--with-tsil-incdir=

--help
--version
"

    # handle --with-addons=
    if [[ ${prev}${cur} == "--with-addons=" || ${pprev}${prev} == "--with-addons=" ]] ; then
        local available_addons=$(__find_available_addons)
        COMPREPLY=( $(__select_from_comma_list "${pprev}" "${prev}" "${cur}" "--with-addons=" "$available_addons") )
        return 0
    fi

    # handle --with-solvers=
    if [[ ${prev}${cur} == "--with-solvers=" || ${pprev}${prev} == "--with-solvers=" ]] ; then
        COMPREPLY=( $(__select_from_comma_list "${pprev}" "${prev}" "${cur}" "--with-solvers=" "all two_scale lattice semi_analytic shooting") )
        return 0
    fi

    # handle --with-loop-libraries=
    if [[ ${prev}${cur} == "--with-loop-libraries=" || ${pprev}${prev} == "--with-loop-libraries=" ]] ; then
        COMPREPLY=( $(__select_from_comma_list "${pprev}" "${prev}" "${cur}" "--with-loop-libraries=" "collier looptools fflite") )
        return 0
    fi

    # handle --with-cppflags=
    if [[ ${prev}${cur} == "--with-cppflags=" || ${pprev}${prev} == "--with-cppflags=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-llflags=
    if [[ ${prev}${cur} == "--with-llflags=" || ${pprev}${prev} == "--with-llflags=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-lllibs=
    if [[ ${prev}${cur} == "--with-lllibs=" || ${pprev}${prev} == "--with-lllibs=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-cxx=
    if [[ ${prev}${cur} == "--with-cxx=" || ${pprev}${prev} == "--with-cxx=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-cxx=" $(compgen -c)) )
        return 0
    fi

    # handle --with-cxxflags=
    if [[ ${prev}${cur} == "--with-cxxflags=" || ${pprev}${prev} == "--with-cxxflags=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-fc=
    if [[ ${prev}${cur} == "--with-fc=" || ${pprev}${prev} == "--with-fc=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-fc=" $(compgen -c)) )
        return 0
    fi

    # handle --with-fflags=
    if [[ ${prev}${cur} == "--with-fflags=" || ${pprev}${prev} == "--with-fflags=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-gsl-config=
    if [[ ${prev}${cur} == "--with-gsl-config=" || ${pprev}${prev} == "--with-gsl-config=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-gsl-config=" $(compgen -c)) )
        return 0
    fi

    # handle --with-math-cmd=
    if [[ ${prev}${cur} == "--with-math-cmd=" || ${pprev}${prev} == "--with-math-cmd=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-math-cmd=" $(compgen -c)) )
        return 0
    fi

    # handle --with-models=
    if [[ ${prev}${cur} == "--with-models=" || ${pprev}${prev} == "--with-models=" ]] ; then
        local available_models=$(__find_created_models)
        COMPREPLY=( $(__select_from_comma_list "${pprev}" "${prev}" "${cur}" "--with-models=" "all $available_models") )
        return 0
    fi

    # handle --with-optional-modules=
    if [[ ${prev}${cur} == "--with-optional-modules=" || ${pprev}${prev} == "--with-optional-modules=" ]] ; then
        COMPREPLY=( $(__select_from_comma_list "${pprev}" "${prev}" "${cur}" "--with-optional-modules=" "examples test") )
        return 0
    fi

    # handle --with-shared-ldflags
    if [[ ${prev}${cur} == "--with-shared-ldflags=" || ${pprev}${prev} == "--with-shared-ldflags=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-shared-lib-ext=
    if [[ ${prev}${cur} == "--with-shared-lib-ext=" || ${pprev}${prev} == "--with-shared-lib-ext=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-shared-lib-cmd=
    if [[ ${prev}${cur} == "--with-shared-lib-cmd=" || ${pprev}${prev} == "--with-shared-lib-cmd=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-shared-lib-cmd=" $(compgen -c)) )
        return 0
    fi

    # handle --with-static-ldflags
    if [[ ${prev}${cur} == "--with-static-ldflags=" || ${pprev}${prev} == "--with-static-ldflags=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-static-lib-ext=
    if [[ ${prev}${cur} == "--with-static-lib-ext=" || ${pprev}${prev} == "--with-static-lib-ext=" ]] ; then
        COMPREPLY=()
        return 0
    fi

    # handle --with-static-lib-cmd=
    if [[ ${prev}${cur} == "--with-static-lib-cmd=" || ${pprev}${prev} == "--with-static-lib-cmd=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-static-lib-cmd=" $(compgen -c)) )
        return 0
    fi

    __build_filename_completion_list "${opts}"
}

_createaddon()
{
    local cur opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    opts="
--name=
--force
--help
"

    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
    if [[ ${#COMPREPLY[@]} == 1 && ${COMPREPLY[0]} != "--"*"=" ]] ; then
        # if there's only one option, without =, then allow a space
        compopt +o nospace
    fi

    return 0
}

_createmodel()
{
    local cur prev pprev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    pprev="${COMP_WORDS[COMP_CWORD-2]}"
    opts="
--name=
--model-file=
--sarah-model=
--with-math-cmd=
--enable-debug
--disable-debug
--force
--help
"

    # handle --model-file=
    if [[ ${prev}${cur} == "--model-file=" || ${pprev}${prev} == "--model-file=" ]] ; then
        local available_models=$(__find_available_models)
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--model-file=" "$available_models") )
        return 0
    fi

    # handle --with-math-cmd=
    if [[ ${prev}${cur} == "--with-math-cmd=" || ${pprev}${prev} == "--with-math-cmd=" ]] ; then
        COMPREPLY=( $(__select_from_list "${pprev}" "${prev}" "${cur}" "--with-math-cmd=" $(compgen -c)) )
        return 0
    fi

    # handle other options
    COMPREPLY=( $(compgen -W "${opts}" -- "${cur}") )
    if [[ ${#COMPREPLY[@]} == 1 && ${COMPREPLY[0]} != "--"*"=" ]] ; then
        # if there's only one option, without =, then allow a space
        compopt +o nospace
    fi

    return 0
}

echo "installing bash completion for configure"
complete -o nospace -F _configure configure

echo "installing bash completion for createmodel"
complete -o nospace -F _createmodel createmodel

echo "installing bash completion for createaddon"
complete -o nospace -F _createaddon createaddon

spectrum_generators=$(find models/ -name "run_*.x" -a ! -name "run_cmd_line*.x" -executable)

for sg in ${spectrum_generators}
do
    echo "installing bash completion for ${sg}"
    complete -o nospace -F _run_spectrum_generator "${sg}"
    complete -o nospace -F _run_spectrum_generator "./${sg}"
done
