digital_header = digital_data.h bit_operator.h powtwo.h config.h

lib_LIBRARIES = libtrmvnorm.a

libtrmvnorm_a_SOURCES = TruncatedNormal.cpp \
	MCQMCIntegration.cpp \
	DigitalNet.cpp $(digital_header) \
	BotevAlgorithm.cpp kahan.hpp \
	pnorm.cpp pnorm.h \
	nxlw_data.cpp \
	solw_data.cpp

AM_CXXFLAGS = -I../include -O3 -Wall -Wextra -D__STDC_CONSTANT_MACROS
