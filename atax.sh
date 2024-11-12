#!/bin/bash

# Компилируем программу
gcc -o atax atax.c -fopenmp -lm -std=c99

# Массив с количеством потоков для тестирования
thread_counts=(1 2 4)

# Запуск программы с каждым количеством потоков
for threads in "${thread_counts[@]}"; do
    echo "Running with $threads threads:"
    # Устанавливаем количество потоков
    export OMP_NUM_THREADS=$threads
    # Запускаем программу и измеряем время выполнения
    ./atax
    echo ""
done
