SHELL=/bin/bash

CPPFLAGS += -std=c++11 -W -Wall  -g
CPPFLAGS += -O3 -lrt
#-ltbb
CPPFLAGS += -I include

CLIENT=src/bitecoin_miner
EXCHANGE_ADDR = 155.198.117.237
EXCHANGE_PORT = 4123
# For your makefile, add TBB and OpenCL as appropriate

# Launch client and server connected by pipes
launch_pipes: src/bitecoin_server $(CLIENT)
	-rm .fifo_rev
	mkfifo .fifo_rev
	# One direction via pipe, other via fifo
	$(CLIENT) client1 3 file .fifo_rev - |\
		(src/bitecoin_server server1 3 file - .fifo_rev &> /dev/null)

# Launch an "infinite" server, that will always relaunch
launch_infinite_server: src/bitecoin_server
	while [ 1 ]; do \
		src/bitecoin_server server1-$USER 3 tcp-server 4000; \
	done;

# Launch a client connected to a local server
connect_local: $(CLIENT)
	$(CLIENT) client-$USER 3 tcp-client localhost 4000

# Launch a client connected to a shared exchange
connect_exchange: $(CLIENT)
	$(CLIENT) client-$(USER) 3 tcp-client $(EXCHANGE_ADDR)  $(EXCHANGE_PORT)

clean:
	rm -f src/*.o src/bitecoin_client src/bitecoin_miner src/bitecoin_server
