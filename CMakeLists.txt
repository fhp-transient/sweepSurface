cmake_minimum_required (VERSION 3.16)
project ("project")

file(GLOB SRCS "./src/*.cpp" "./src/*.c" "./imgui/*.cpp" "./imgui/backends/*.cpp")
add_executable(project ${SRCS}
        include/Algorithm.h
        include/myutil.h
)
target_include_directories(project PRIVATE "./include" "./imgui" "./imgui/backends")
if(LINUX)
    target_link_libraries(project dl GL glfw)
else()
    target_link_directories(project PRIVATE "./lib")
    target_link_libraries(project glfw3)
endif()
