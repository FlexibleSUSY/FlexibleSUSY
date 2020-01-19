#! /bin/sh

# This script returns the include paths for Mathematica's MathLink and
# LibraryLink .
# Author: Alexander Voigt

math_cmd=math

find_math_dirs() {
    eval eval `printf "%s" '
       Print["sysid=\"", $SystemID, "\""];
       Print["topdir=\"", $TopDirectory, "\""];
       Exit[]' | "${math_cmd}" -noprompt | tr '\r' ' ' | tail -2`

    # check whether Cygwin's dlltool can handle 64-bit DLLs
    test "$sysid" = Windows-x86-64 && {
        ${DLLTOOL:-dlltool} --help | grep x86-64 > /dev/null || sysid=Windows
    }

    topdir=`cd "$topdir" ; echo $PWD`
}

get_librarylink_incpath() {
    find_math_dirs

    for p in \
        "$topdir/SystemFiles/IncludeFiles/C" ; do
        test -d "$p" && break
    done

    echo "$p"
}

get_mathlink_incpath() {
    find_math_dirs

    for p in \
        "$topdir/SystemFiles/Links/MathLink/DeveloperKit/$sysid/CompilerAdditions" \
        "$topdir/SystemFiles/Links/MathLink/DeveloperKit/CompilerAdditions" \
        "$topdir/AddOns/MathLink/DeveloperKit/$sysid/CompilerAdditions" ; do
        test -d "$p" && break
    done

    echo "$p"
}

usage() {
cat <<EOF
Usage: ./$(basename $0) [options]
Options:

  -I                Prepend -I on each following include path.
  --librarylink,-l  Output include path for LibraryLink headers.
  --mathlink,-m     Output include path for MathLink headers.
  --math-cmd=       Use Mathematica kernel.
  --help,-h         Print this help message.
EOF
}

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=$(echo "$1" | sed 's/[-_a-zA-Z0-9]*=//') ;;
            *) optarg= ;;
        esac

        case $1 in
            -I)               prepend="-I" ;;
            --librarylink|-l) path="${path} '${prepend}$(get_librarylink_incpath)'" ;;
            --math-cmd=*)     math_cmd=$optarg ;;
            --mathlink|-m)    path="${path} '${prepend}$(get_mathlink_incpath)'" ;;
            --help|-h)        usage; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

echo "$path"
