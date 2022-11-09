CC = g++
CFLAGS= -std=c++11 -Wall -pedantic -Werror
BIN_DIR = bin

ES1_1 = $(BIN_DIR)/Errori_1
ES1_2 = $(BIN_DIR)/Errori_2
ES1_3 = $(BIN_DIR)/Errori_3
ES2_1 = $(BIN_DIR)/Sistemi_1
ES2_2 = $(BIN_DIR)/Sistemi_2
ES2_3 = $(BIN_DIR)/Sistemi_3

ES1_1_OBJS = $(BIN_DIR)/Errori_1.o 
ES1_2_OBJS = $(BIN_DIR)/Errori_2.o 
ES1_3_OBJS = $(BIN_DIR)/Errori_3.o 
ES2_1_OBJS = $(BIN_DIR)/Sistemi_1.o 
ES2_2_OBJS = $(BIN_DIR)/Sistemi_2.o 
ES2_3_OBJS = $(BIN_DIR)/Sistemi_3.o 

EXECS = $(ES1_1) $(ES1_2) $(ES1_3) $(ES2_1) $(ES2_2) $(ES2_3)

all: $(EXECS)

.PHONY: clean tgz

# Esercizi Errori
$(ES1_1): $(ES1_1_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

$(ES1_2): $(ES1_2_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

$(ES1_3): $(ES1_3_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

$(BIN_DIR)/Errori_1.o: Errori/Es1-1.cpp | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ -c Errori/Es1-1.cpp

$(BIN_DIR)/Errori_2.o: Errori/Es1-2.cpp | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ -c Errori/Es1-2.cpp

$(BIN_DIR)/Errori_3.o: Errori/Es1-3.cpp | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ -c Errori/Es1-3.cpp



# Esercizi Sistemi Lineari
$(ES2_1): $(ES2_1_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

$(ES2_2): $(ES2_2_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

$(ES2_3): $(ES2_3_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

$(BIN_DIR)/Sistemi_1.o: Sistemi/Es2-1.cpp | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ -c Sistemi/Es2-1.cpp

$(BIN_DIR)/Sistemi_2.o: Sistemi/Es2-2.cpp | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ -c Sistemi/Es2-2.cpp

$(BIN_DIR)/Sistemi_3.o: Sistemi/Es2-3.cpp | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ -c Sistemi/Es2-3.cpp


# Utilities
clean:
	rm -rf $(EXECS) $(BIN_DIR)/*.o

tgz: clean
		cd ..; tar cvzf Errori_Sistemi.tgz Errori_Sistemi
