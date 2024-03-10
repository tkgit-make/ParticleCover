CXX=g++
RM=rm -f
CPPFLAGS=-g
LDFLAGS=-g
LDLIBS=
SRCS=makePatch_alignedtoLine.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

makePatch: $(OBJS)
	 $(CXX)  $(LDFLAGS) -o makePatch $(OBJS) $(LDLIBS)

makePatch_alignedtoLine.o: makePatch_alignedtoLine.cpp
	g++ $(CPPFLAGS) -c makePatch_alignedtoLine.cpp

clean: 
	$(RM) $(OBJS)
