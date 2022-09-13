#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

config_dir="${BASEDIR}/../config"
space_dir="${BASEDIR}/space dir"
space_dir_touched_file="${space_dir}/touched"
space_dir_makefile="${space_dir}/Makefile"
space_dir_makefile_tmpl="${space_dir}/Makefile.in"
error=0

cleanup() {
    [ -e "${space_dir_makefile}" ] && (cd "$space_dir" && make clean)
    rm -f "${space_dir_touched_file}" "${space_dir_makefile}"
}

create_makefile() {
    sed -e "s|@CONFIG_DIR@|${config_dir}|" \
        < "${space_dir_makefile_tmpl}" > "${space_dir_makefile}"
}

cleanup
create_makefile

(cd "${space_dir}" && make)

if [ -e "${space_dir_touched_file}" ]; then
    echo "Info: file found: ${space_dir_touched_file}"
else
    echo "Error: file not found: ${space_dir_touched_file}"
    error=1 
fi

if [ -e "${space_dir_touched_file}_1" ]; then
    echo "Info: file found: ${space_dir_touched_file}_1"
else
    echo "Error: file not found: ${space_dir_touched_file}_1"
    error=1
fi
if [ -e "${space_dir_touched_file}_2" ]; then
    echo "Info: file found: ${space_dir_touched_file}_2"
else
    echo "Error: file not found: ${space_dir_touched_file}_2"
    error=1
fi

cleanup

exit ${error}
