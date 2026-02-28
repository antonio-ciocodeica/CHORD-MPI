#!/bin/bash

echo "==========================================="
echo "   Se ruleaza toate testele de CHORD MPI   "
echo "==========================================="

CHECKER=check_lookup.py
TIMEOUT=10

### Scoruri inițiale ###
CORRECT_SCORE=0
EFF_SCORE=0

###  Pasul de build
echo ""
echo ">>> Se compilează proiectul cu make build ..."
pushd ../src > /dev/null || exit 1
make build > build_log.txt 2>&1
BUILD_STATUS=$?

if [ $BUILD_STATUS -ne 0 ]; then
    echo "Eroare: compilarea a eșuat!"
    cat build_log.txt
    echo ""
    echo "Corectitudine: 0/40"
    echo "Eficiență:     0/40"
    echo "Total:         0/80"
    exit 0
fi

echo "Compilare reușită."
echo ""

rm -rf build_log.txt
popd > /dev/null || exit 1

cp ../src/chord .

### Se detectează tool-ul de timeout ###
if command -v timeout >/dev/null 2>&1; then
    TIMEOUT_CMD=timeout
elif command -v gtimeout >/dev/null 2>&1; then
    TIMEOUT_CMD=gtimeout
elif command -v /usr/local/opt/coreutils/libexec/gnubin/timeout >/dev/null 2>&1; then
    TIMEOUT_CMD=/usr/local/opt/coreutils/libexec/gnubin/timeout
else
    echo "Eroare: nu există 'timeout' sau 'gtimeout'."
    exit 1
fi

### Funcție care rulează un test și calculează scorul ###
run_test() {
    TESTDIR=$1
    echo ""
    echo "--- Se rulează $TESTDIR ---"

    NP=$(ls $TESTDIR/in*.txt | wc -l)
    rm -f in*.txt out.txt
    cp $TESTDIR/in*.txt .

    # rulare cu timeout
    $TIMEOUT_CMD $TIMEOUT mpirun --oversubscribe --allow-run-as-root -np $NP ./chord > out.txt 2> run_error.log
    EXIT_CODE=$?

    # dacă e timeout
    if [ $EXIT_CODE -eq 124 ]; then
        echo "FAIL: Programul s-a blocat (timeout)"
        rm -f in*.txt out.txt run_error.log
        return
    fi

    FAIL=0
    LOCAL_CORRECT=1
    LOCAL_EFF=1

    while read -r line; do
        key=$(echo $line | awk '{print $1}')
        succ=$(echo $line | awk '{print $2}')
        max_hops=$(echo $line | awk '{print $3}')

        out_line=$(grep "Lookup $key:" out.txt)

        if [ -z "$out_line" ]; then
            echo "FAIL: Lipsește lookup pentru cheia $key"
            LOCAL_CORRECT=0
            LOCAL_EFF=0
            continue
        fi

        # se rulează checker-ul
        python3 $CHECKER "$line" "$out_line"
        RESULT=$?

        if [ $RESULT -eq 2 ]; then
            # succesor greșit
            LOCAL_CORRECT=0
            LOCAL_EFF=0
        elif [ $RESULT -eq 3 ]; then
            # prea multe hopuri
            LOCAL_EFF=0
        fi

    done < $TESTDIR/expected.txt

    # update scoruri
    if [ $LOCAL_CORRECT -eq 1 ]; then
        CORRECT_SCORE=$((CORRECT_SCORE + 5))
    fi
    if [ $LOCAL_EFF -eq 1 ]; then
        EFF_SCORE=$((EFF_SCORE + 5))
    fi

    # afișare output test
    if [ $LOCAL_CORRECT -eq 1 ] && [ $LOCAL_EFF -eq 1 ]; then
        echo "--- $TESTDIR: PASS TOTAL (10/10) ---"
    elif [ $LOCAL_CORRECT -eq 1 ] && [ $LOCAL_EFF -eq 0 ]; then
        echo "--- $TESTDIR: PASS PARTIAL (5/10) ---"
    else
        echo "--- $TESTDIR: FAIL (0/10) ---"
    fi

    rm -f in*.txt out.txt run_error.log
}

### Se rulează testele ###
run_test tests/test1
run_test tests/test2
run_test tests/test3
run_test tests/test4
run_test tests/test5
run_test tests/test6
run_test tests/test7
run_test tests/test8

make clean > /dev/null 2>&1

### Rezultate finale ###
echo ""
echo "--------------------------------------"
echo ""
echo "Corectitudine: $CORRECT_SCORE/40"
echo "Eficiență:     $EFF_SCORE/40"
echo "Total:         $((CORRECT_SCORE + EFF_SCORE))/80"
