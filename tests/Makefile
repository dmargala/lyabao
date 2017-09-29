test: test_trig1.txt compare_distance.txt test_trig_examples

test_trig_examples: test_trig_examples.cpp
	g++ test_trig_examples.cpp -lcosmo -o test_trig_examples

test_trig1.txt: test_trig.py
	python test_trig.py

compare_distance.txt: test_cosmo_distance.py test_cosmo_distance.txt
	python test_cosmo_distance.py

test_cosmo_distance.txt: test_cosmo_distance
	./test_cosmo_distance > test_cosmo_distance.txt

test_cosmo_distance: test_cosmo_distance.cpp
	g++ test_cosmo_distance.cpp -lcosmo -o test_cosmo_distance

clean:
	rm test_trig?.txt \
		compare_distance.txt  \
		test_cosmo_distance.txt \
		test_cosmo_distance \
		test_trig_examples


.PHONY: clean
