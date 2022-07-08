all: build/main.o build/complex.o build/vector_arm64.o
	g++ build/main.o build/complex.o build/vector_arm64.o -o fft -O3

build/main.o: src/main.cpp headers/complex.hpp
	g++ -c src/main.cpp -o build/main.o -O3

build/vector_arm64.o: src/vector_arm64.cpp headers/vector_ops.hpp 
	g++ -c src/vector_arm64.cpp -o build/vector_arm64.o -O3

build/complex.o: src/complex.cpp headers/complex.hpp headers/vector_ops.hpp
	g++ -c src/complex.cpp -o build/complex.o -O3
clean:
	rm build/*.o fft