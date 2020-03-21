include_directories(${PROJECT_SOURCE_DIR}/ext)
message(${PROJECT_SOURCE_DIR}/ext)

file(GLOB src_file
        *.h
        *.cpp)

add_library( lemath ${src_file})