CXX=g++-11
RM=rm -f
CPPFLAGS=-g -std=c++20
LDFLAGS=-g
LDLIBS=
SRCS=makePatch_alignedtoLine.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

makePatch: $(OBJS)
	 $(CXX)  $(LDFLAGS) -o makePatch $(OBJS) $(LDLIBS)

makePatch_alignedtoLine.o: makePatch_alignedtoLine.cpp
	$(CXX) $(CPPFLAGS) -c makePatch_alignedtoLine.cpp

clean: 
	$(RM) $(OBJS)
