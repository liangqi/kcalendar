project (test)
 
find_package(KDE4 REQUIRED)
include_directories( ${KDE4_INCLUDES} )
 
set(test_SRCS main.cpp)
 
kde4_add_executable(test ${test_SRCS})
target_link_libraries(test ${KDE4_KDEUI_LIBS})

