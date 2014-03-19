SHELL=/bin/bash

CPPFLAGS += -std=c++11 -W -Wall  -g
CPPFLAGS += -O3
CPPFLAGS += -I include -I src/miner
LDLIBS += -lrt -ltbb

CLIENT=src/bitecoin_miner
EXCHANGE_ADDR = 155.198.117.237
EXCHANGE_PORT = 4123
LOG_LEVEL ?= 3
CLIENT_ID = $(shell head -c 512 /dev/urandom | md5sum | cut -c 1-10)
#CLIENT_ID="Donkey++"

# For your makefile, add TBB and OpenCL as appropriate

all: src/bitecoin_server src/bitecoin_miner src/bitecoin_client

# Launch client and server connected by pipes
launch_pipes: src/bitecoin_server $(CLIENT)
	-rm .fifo_rev
	mkfifo .fifo_rev
	# One direction via pipe, other via fifo
	$(CLIENT) client1 $(LOG_LEVEL) file .fifo_rev - |\
		(src/bitecoin_server server1 $(LOG_LEVEL) file - .fifo_rev &> /dev/null)

# Launch an "infinite" server, that will always relaunch
launch_infinite_server: src/bitecoin_server
	while [ 1 ]; do \
		src/bitecoin_server server1-$USER $(LOG_LEVEL) tcp-server 4000; \
	done;

launch_server: src/bitecoin_server
	src/bitecoin_server server1-$USER $(LOG_LEVEL) tcp-server 4000;

# Launch a client connected to a local server
connect_local: $(CLIENT)
	$(CLIENT) client-$USER $(LOG_LEVEL) tcp-client localhost 4000

# Launch a client connected to a shared exchange
connect_exchange: $(CLIENT)
	$(CLIENT) $(CLIENT_ID) $(LOG_LEVEL) tcp-client $(EXCHANGE_ADDR)  $(EXCHANGE_PORT)

clean:
	rm -f src/*.o src/bitecoin_client src/bitecoin_miner src/bitecoin_server\
		src/miner/*.o

.PHONY: launch_pipes launch_infinite_server connect_exchange connect_local\
	clean all launch_server
