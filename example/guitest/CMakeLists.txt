project (guitest)
 
find_package(KDE4 REQUIRED)
include_directories( ${KDE4_INCLUDES} )
 
set(guitest_SRCS main.cpp)
 
kde4_add_executable(guitest ${guitest_SRCS})
target_link_libraries(guitest ${KDE4_KDEUI_LIBS})

