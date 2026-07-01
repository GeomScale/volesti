# Spectrahedron / SDPA format tests

add_executable(sdpa_format_test sdpa_format_test.cpp $<TARGET_OBJECTS:test_main>)
add_test(NAME test_sdpa_raw_roundtrip COMMAND sdpa_format_test -tc="sdpa_raw_roundtrip")
add_test(NAME test_sdpa_spectrahedron_roundtrip COMMAND sdpa_format_test -tc="sdpa_spectrahedron_roundtrip")
add_test(NAME test_sdpa_comment_lines COMMAND sdpa_format_test -tc="sdpa_comment_lines")
add_test(NAME test_sdpa_diagonal_blocks COMMAND sdpa_format_test -tc="sdpa_diagonal_blocks")
add_test(NAME test_sdpa_error_truncated COMMAND sdpa_format_test -tc="sdpa_error_truncated_file")
TARGET_LINK_LIBRARIES(sdpa_format_test coverage_config)

# test_sdpa_format.cpp is stale — uses deprecated API, not built
