add_executable(lpmd lpmd_calc.cpp)
target_link_libraries(lpmd Eigen3::Eigen)
install(TARGETS lpmd DESTINATION bin)
