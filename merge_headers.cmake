message(STATUS "Call merge_headers")
message(STATUS "Writing to ${OUTPUT}")
message(STATUS "Merging from ${MERGE_LIST}")
message(STATUS "Prefix : '${PREFIX}'")

file(WRITE ${OUTPUT} "#pragma once\n\n")
file(STRINGS ${MERGE_LIST} cur_dir_files)
foreach(SUBFILE ${cur_dir_files})
	get_filename_component(FILENAME ${SUBFILE} NAME)
	file(APPEND ${OUTPUT} "#include \"" ${PREFIX}${FILENAME} "\"\n")
endforeach()
