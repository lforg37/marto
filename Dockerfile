FROM alpine:3.10.3
RUN apk add bash git g++ cmake make boost-dev mpfr-dev gmp-dev; \
	git clone https://github.com/yuguen/hint.git ; \
	cd hint && mkdir build && cd build && cmake .. && make install && cd ../..; \
	git clone https://github.com/Xilinx/HLS_arbitrary_Precision_Types.git ap_int; \

	git clone https://gitlab.com/cerlane/SoftPosit.git ; \
	cd SoftPosit/build/Linux-x86_64-GCC && make -j6 all ; \

	cd ../../..; \
	git clone https://gitlab.inria.fr/lforget/marto.git ; \
	cd marto && mkdir build && cd build ; \
	cmake .. -DSOFTPOSIT_H=/SoftPosit/source/include \
         -DSOFTPOSIT_LIB=/SoftPosit/build/Linux-x86_64-GCC/softposit.a \
         -DAP_INT_INCLUDE_DIR=/ap_int/include \
         -DBUILD_UNIT_TEST=1 \
         -DCMAKE_BUILD_TYPE=RELEASE ; \
    make -j ; \


    ./testLibNumForm_exe ; \
    ./testIEEEAdderGMP_exe ; \
    ./testIEEEMultiplier_exe ; \
    ./testAritOps_exe -t @long ; \
    ./testConvert_exe ; \
    ./testKulisch_exe -t @long ; \