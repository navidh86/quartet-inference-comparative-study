echo "" > progress_$1.txt
for (( i=$1; i<=$2; i++ )); do
    echo $i
    python3 SCRIPT_RUNNER_Avian.py $i
    echo $i >> progress_$1.txt
done
echo "done##########################################"

# for i in {5..9}; do
#     echo $i
#     python3 SCRIPT_RUNNER_Avian.py $i &
# done