

walks: graphics_for_walks.o 
	g++ $(CFLAGS) -o walks graphics_for_walks.o -lGL -lm -lgslcblas -lgsl -lgslcblas -lgsl -lglut

graphics_for_walks.o: graphics_for_walks.cpp
	g++ $(CFLAGS) -c graphics_for_walks.cpp

clean:
	rm -f *.o walks
