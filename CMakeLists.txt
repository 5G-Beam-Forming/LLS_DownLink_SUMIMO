cmake_minimum_required (VERSION 2.8.11)
project (DL_SUMIMO.exe)

add_executable (DL_SUMIMO
src/LDPC/method_ChannelCoding.cpp
src/LDPC/RX_AWGN_Channel.cpp
src/LDPC/Rx_ChannelDecoding.cpp
src/LDPC/TBS_determination.cpp
src/LDPC/Tx_ChannelEncoding.cpp
src/LDPC/method_ChannelCoding.cpp
src/LDPC/method_ChannelDecoding.cpp
src/LDPC/method_ChannelEncoding.cpp
src/LDPC/method_QAM_Modulation.cpp
src/LDPC/misc.cpp
src/ChannelOutput_CDL.cpp
src/ChannelOutput_TDL.cpp
src/getCQITable.cpp
src/LDPC_CC.CPP
src/main_SUMIMO.cpp
src/module_BS_SUMIMO.cpp
src/module_Channel_5G.cpp
src/module_Parameter_MIMO.cpp
src/module_UE_MIMO.cpp
src/sumimo.cpp)

target_link_libraries (DL_SUMIMO LINK_PUBLIC armadillo)