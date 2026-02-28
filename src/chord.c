#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define M 4 // 2^m = 16 ID-uri posibile
#define RING_SIZE 16
#define MAX_NODES 32
#define MAX_PATH 32

/************************************************************
 * Structuri CHORD
 ************************************************************/
typedef struct {
    int start;              // start[i] = (id + 2^i) mod 2^m
    int node;               // finger[i] = succesor(start)
} Finger;

typedef struct {
    int id;                 // ID CHORD al nodului curent
    int successor;          // succesorul în inel
    int predecessor;        // predecesorul în inel
    Finger finger[M];       // finger table-ul (de completat de voi)
} NodeState;

typedef struct {
    int initiator_id;       // cine a inițiat lookup-ul
    int current_id;         // nodul care procesează mesajul
    int seq;                // numarul de secventa al mesajului (folosit pentru a pastra ordinea in afisare)
    int key;                // cheia căutată
    int path[MAX_PATH];     // traseul lookup-ului
    int path_len;           // lungimea traseului
} LookupMsg;

/************************************************************
 * Tag-uri mesaje MPI
 ************************************************************/
enum {
    TAG_LOOKUP_REQ = 1,
    TAG_LOOKUP_REP,
    TAG_DONE
};

/************************************************************
 * Variabile globale utile
 ************************************************************/
int world_rank, world_size;
NodeState self;

int all_ids[MAX_NODES];      // ID CHORD pentru fiecare rank
int id_to_rank[RING_SIZE];   // mapare inversă: CHORD ID -> MPI rank
int sorted_ids[MAX_NODES];   // ID-uri sortate ale nodurilor din inel
int num_nodes;

/************************************************************
 * Funcții utile pentru interval circular CHORD
 ************************************************************/
int in_interval(int x, int start, int end) {
    if (start < end)
        return (x > start && x <= end);
    if (start > end)
        return (x > start || x <= end);
    return 1;   // intervalul acoperă tot cercul
}

/************************************************************
 * Construirea mapării rank -> id și id -> rank
 ************************************************************/
void build_id_maps() {
    for (int i = 0; i < RING_SIZE; i++)
        id_to_rank[i] = -1;

    for (int r = 0; r < world_size; r++) {
        int id = all_ids[r];
        if (id >= 0 && id < RING_SIZE)
            id_to_rank[id] = r;
    }
}

int rank_from_id(int id) {
    if (id < 0 || id >= RING_SIZE)
        return -1;
    return id_to_rank[id];
}

/************************************************************
 * Construirea inelului CHORD static
 ************************************************************/
int cmp_int(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

int cmp_msg(const void *a, const void *b) {
    LookupMsg* msg_a = (LookupMsg*) a; 
    LookupMsg* msg_b = (LookupMsg*) b;
    return msg_a->seq - msg_b->seq;
}

void build_global_ring() {
    num_nodes = world_size;

    for (int i = 0; i < num_nodes; i++)
        sorted_ids[i] = all_ids[i];

    qsort(sorted_ids, num_nodes, sizeof(int), cmp_int);

    for (int i = 0; i < num_nodes; i++) {
        int id = sorted_ids[i];
        int succ = sorted_ids[(i + 1) % num_nodes];
        int pred = sorted_ids[(i - 1 + num_nodes) % num_nodes];

        if (id == self.id) {
            self.successor = succ;
            self.predecessor = pred;
        }
    }
}

/************************************************************
 * Construirea finger table-ului.
 * find_successor_simple() caută succesorul lui start[i] în sorted_ids[],
 * care conține doar nodurile existente.
 ************************************************************/
int find_successor_simple(int key) {
    for (int i = 0; i < num_nodes; i++)
        if (sorted_ids[i] >= key)
            return sorted_ids[i];
    return sorted_ids[0]; // wrap-around
}

/******************************************************
 * Succesorul se caută doar în lista sorted_ids.
 * Inelul este static, deci finger table-ul este static.
 ******************************************************/
void build_finger_table() {
    for (int i = 0; i < M; i++) {
        int pow = 1 << i;
        int start = (self.id + pow) % RING_SIZE;
        int node = find_successor_simple(start);
        
        self.finger[i].start = start;
        self.finger[i].node = node;
    }
}

/************************************************************
 *  Această funcție returnează cel mai "mare" finger care
 *  se află strict în intervalul CHORD (self.id, key).
 *  Dacă niciun finger nu se potrivește, se intoarce succesorul.
 ************************************************************/
int closest_preceding_finger(int key) {
    for (int i = M - 1; i >= 0; i--) {
        int finger = self.finger[i].node;

        // Interval deschis
        if (self.id == finger || finger == key) {
            continue;
        }

        if (in_interval(finger, self.id, key)) {
            return finger;
        }
    }

    return self.successor;   // fallback
}

/************************************************************
 *   Aceasta este funcția esențială pentru rutare distribuită.
 *   Pașii simplificati CHORD :
 *
 *   1. Se adauga self.id în path.
 *   2. Dacă succesorul nostru este responsabil de key:
 *           - se adauga succesorul în path
 *           - se trimite TAG_LOOKUP_REP către inițiator
 *   3. Altfel:
 *           - next = closest_preceding_finger(key)
 *           - se trimite TAG_LOOKUP_REQ către closest_preceding_finger
 *
 ************************************************************/
void handle_lookup_request(LookupMsg *msg) {
    msg->path[msg->path_len++] = self.id;
    
    if (in_interval(msg->key, self.id, self.successor)) {
        msg->path[msg->path_len++] = self.successor;
        int initiator_rank = rank_from_id(msg->initiator_id);
        MPI_Send(msg, sizeof(LookupMsg), MPI_BYTE, initiator_rank, TAG_LOOKUP_REP, MPI_COMM_WORLD);
    } else {
        int next_node = closest_preceding_finger(msg->key);
        int next_rank = rank_from_id(next_node);
        MPI_Send(msg, sizeof(LookupMsg), MPI_BYTE, next_rank, TAG_LOOKUP_REQ, MPI_COMM_WORLD);
    }
}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /******************************************************
     * Citire input
     ******************************************************/
    char fname[32];
    sprintf(fname, "in%d.txt", world_rank);
    FILE *f = fopen(fname, "r");
    if (!f) {
        fprintf(stderr, "Rank %d: cannot open %s\n", world_rank, fname);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int nr_lookups;
    fscanf(f, "%d", &self.id);
    fscanf(f, "%d", &nr_lookups);

    int *lookups = malloc(nr_lookups * sizeof(int));
    for (int i = 0; i < nr_lookups; i++)
        fscanf(f, "%d", &lookups[i]);
    fclose(f);

    /******************************************************
     * Distribuirea ID-urilor tuturor nodurilor
     ******************************************************/
    MPI_Allgather(&self.id, 1, MPI_INT,
                  all_ids, 1, MPI_INT,
                  MPI_COMM_WORLD);

    build_id_maps();
    build_global_ring();
    build_finger_table();

    MPI_Barrier(MPI_COMM_WORLD);

    /************************************************************
     * Inițierea lookup-urilor locale
     *
     *    Pentru fiecare valoare Key citită din input:
     *       - se construieste LookupMsg
     *       - se trimite un mesaj de tip TAG_LOOKUP_REQ către propriul rank
     ************************************************************/

    for (int i = 0; i < nr_lookups; i++) {
        LookupMsg* msg = malloc(sizeof(LookupMsg));

        msg->initiator_id = self.id;
        msg->current_id = self.id;
        msg->seq = i;
        msg->key = lookups[i];
        msg->path_len = 0;

        int rank = rank_from_id(self.id);
        MPI_Send(msg, sizeof(LookupMsg), MPI_BYTE, rank, TAG_LOOKUP_REQ, MPI_COMM_WORLD);
    }

    free(lookups);

    /************************************************************
     * Service-loop distribuit
     *   - se primesc mesaje cu MPI_Recv()
     *   - se proceseaza TAG_LOOKUP_REQ → handle_lookup_request()
     *   - se proceseaza TAG_LOOKUP_REP → afișare + countdown local
     *   - se trimite TAG_DONE tuturor când se termina
     *   - se opreste loop-ul doar când se primeste DONE de la toate nodurile
     ************************************************************/

    LookupMsg* messages = malloc(nr_lookups * sizeof(LookupMsg));
    int lookups_done = 0;
    int done_received = 0;
    int is_done = 0;

    while (1) {
        MPI_Status status;
        LookupMsg msg;

        MPI_Recv(&msg, sizeof(LookupMsg), MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        // Proceseaza requesturile
        if (status.MPI_TAG == TAG_LOOKUP_REQ) {
            handle_lookup_request(&msg);
        }

        // Adauga raspunsurile in vectorul de mesaje pe pozitia data de numarul de secventa,
        // astfel se va pastra ordinea mesajelor citite din fisierele input
        if (status.MPI_TAG == TAG_LOOKUP_REP) {
            messages[msg.seq] = msg;
            lookups_done++;
        }

        // Incrementeaza numarul de mesaje TAG_DONE primite
        if (status.MPI_TAG == TAG_DONE) {
            done_received++;
        }

        // Daca nodul si-a terminat toate lookup-urile, trimite TAG_DONE la toate nodurile
        if (!is_done && lookups_done == nr_lookups) {
            is_done = 1;
            int x;

            for (int i = 0; i < world_size; i++) {
                if (world_rank != i) {
                    MPI_Send(&x, 1, MPI_INT, i, TAG_DONE, MPI_COMM_WORLD);
                }
            }
        }

        // Daca s-a primit TAG_DONE de la toti, termina rutarea
        if (is_done && done_received == world_size - 1) {
            break;
        }
    }

    // Afiseaza path-urile de lookup
    for (int i = 0; i < nr_lookups; i++) {
        LookupMsg msg = messages[i];

        printf("Lookup %d: ", msg.key);
    
        for (int i = 0; i < msg.path_len - 1; i++) {
            printf("%d -> ", msg.path[i]);
        }
    
        printf("%d\n", msg.path[msg.path_len - 1]);
    }

    free(messages);
    MPI_Finalize();
    return 0;
}
