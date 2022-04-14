# функция проверки взаиморасположения праймеров пары
def coord_match(seqs_dict, specific_primer):
    for group in specific_primer:   # для каждой группы
        match = []   # список пар праймеров
        for i in range(len(specific_primer[group])):   # для каждого праймера из списка группы
            primer_1 = specific_primer[group][i]
            for j in range(i + 1, len(specific_primer[group])):   # проверяем каждый следующий праймер этой группы
                primer_2 = specific_primer[group][j]
                for group_num in seqs_dict:   # для каждой группы словаря последовательностей
                    if group_num == group:   # если ключи словарей праймеров и посл-тей совпадают
                        for original_seq in seqs_dict[group_num]:   # для каждой цепи
                            if 'matr' in original_seq.name:   # если оня матричная
                                s_strand = original_seq.seq   # в переменную помещаем ее посл-ть
                                as_strand = original_seq.complement()   # в др переменную помещаем ей компл-ую
                                print(s_strand, as_strand)
                                if primer_1.name == 'forward' and primer_2.name == 'reversed':
                                    primer_2_rev = Primer('reversed', primer_2.reverse())  # переворачиваем обратный пр-р
                                    coord1 = s_strand.index(primer_1.seq)   # находим индекс 1 эл-та прямого праймера
                                    coord2 = as_strand.index(primer_2_rev.seq)   # находим индекс 1 эл-та обратного праймера

                                    print(primer_1.seq, primer_2.seq, primer_2_rev.seq)
                                    if coord1 != coord2:   # если координаты не совпадают, записываем в сет
                                        match.append(primer_1)
                                        match.append(primer_2)
        specific_primer[group] = match   # сохраняем в словаре только отсортированные праймеры

    for group in specific_primer:  # для каждой группы
         for primer in specific_primer[group]:  # для каждого праймера в списке группы
            print(group, primer.name, primer.seq, primer.temp())  # выводим: № группы, тип праймера, температуру отжига

    return specific_primer