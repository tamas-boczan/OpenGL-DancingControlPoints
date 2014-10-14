cmake_minimum_required(VERSION 2.8.4)
project(Grafika)

find_package(GLUT)
find_package(OpenGL)
INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIR})
ADD_EXECUTABLE(grafika Source.cpp)
TARGET_LINK_LIBRARIES(grafika ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pedantic -Wall -Wextra -DDEBUG -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pedantic -Wall -Wextra -g0 -o3")

#+ param, ha nem működne: -Wl,--no-as-needed
# szükséges volt: sudo apt-get install libxmu-dev libxi-dev

# Command line debug: g++ Source.cpp -pedantic -Wall -Wextra -DDEBUG -g -lglut -lGLU -lGL -o grafika
# command line release: g++ Source.cpp -pedantic -Wall -Wextra -g0 -03 -lglut -lGLU -lGL -o grafika