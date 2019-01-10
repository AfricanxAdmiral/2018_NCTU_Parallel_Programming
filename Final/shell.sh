
# makenum
#./makenum
#nvcc fft_cuda.cu -o fft_cuda 
g++ -fopenmp fft.cpp -o fft
g++ -fopenmp fft_serial.cpp -o fft_serial
g++ -fopenmp -lpthread fft_pthread_test.cpp -o fft_pthread_test
g++ mean.cpp -o mean -lm
g++ makenum.cpp -o makenum
g++ makenum_num.cpp -o makenum_num
#./fft s.txt t.txt > fft_output
#./fft_serial s.txt t.txt > fft_serial_output
#diff fft_output fft_serial_output
#rm -f time.txt
#rm -f time_serial.txt
rm -f time_pthread.txt
rm -f time.txt
rm -f time_serial.txt
#rm -f time_cuda.txt

rm -f pthread_size_time.txt
rm -f openmp_size_time.txt
rm -f serial_size_time.txt
#rm -f cuda_size_time.txt
i=1
k=2
while [ "${i}" != "25" ]
do
	echo $i
	./makenum $k $k
	j=0
	while [ "${j}" != 5 ]
	do
		./fft_pthread_test s.txt t.txt 4
		./fft s.txt t.txt 
		./fft_serial s.txt t.txt
		#./fft_cuda s.txt t.txt
		j=$(($j+1))
	done  
	./mean
 	rm -f time_pthread.txt
	rm -f time.txt
	rm -f time_serial.txt
	#rm -f time_cuda.txt
	
	./makenum_num $k
        j=0
	while [ "${j}" != 5 ]
	do
        	./fft s.txt t.txt
		./fft_pthread_test s.txt t.txt
		./fft_serial s.txt t.txt
		#./fft_cuda s.txt t.txt
		j=$(($j+1))
	done
	./mean
	rm -f time_pthread.txt
	rm -f time.txt
	rm -f time_serial.txt
	#rm -f time_cuda.txt
        k=$(($k+1))
        i=$(($i+1))
done

