
for ((k = 1; k <= 32; k *=2)); do
    find -name nkithr.med.$k.* | xargs sed '3~2d' | sort -n -t , -k 1,2 > nkithr.k$k.med.log
    find -name nkithr.min.$k.* | xargs sed '3~2d' | sort -n -t , -k 1,2 > nkithr.k$k.min.log
    find -name nkithr.max.$k.* | xargs sed '3~2d' | sort -n -t , -k 1,2 > nkithr.k$k.max.log
done
