CU_APPS=dudx


all: ${CU_APPS}

%: %.cu
	nvcc -o $@ $<
clean:
	rm -f ${CU_APPS}
