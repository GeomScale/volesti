QD_CPPFLAGS=$(CPPFLAGS) -I$(R_INCLUDE_DIR) -march=native

QD_SOURCES= bits.cc c_dd.cc c_qd.cc dd_const.cc dd_real.cc fpu.cc \
                 qd_const.cc qd_real.cc util.cc

QD_OBJECTS=$(QD_SOURCES:.cc=.o)

libqd.a: $(QD_OBJECTS)
	$(AR) rc libqd.a $(QD_OBJECTS)

.cc.o:
	$(CC) $(CFLAGS) $(CPICFLAGS) $(QD_CPPFLAGS) -c $< -o $@

clean:
	rm -rf $(QD_OBJECTS) libqd.a
