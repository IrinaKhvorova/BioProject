import csv
from Bio import SeqIO
from itertools import zip_longest

class Sequence:

    alphabet = []

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def __getitem__(self, item):  # возвращает символ в последовательности
        return self.seq[item]

    def __len__(self):  # возвращает длину последовательности
        return len(self.seq)

    def name(self):  # возвращает имя последовательности
        return self.name

    @property
    def sequence(self):  # возвращает саму последовательность
        return self.seq


class DNA(Sequence):

    alphabet = ['A', 'T', 'G', 'C']

    def __init__(self, name, seq):
        super().__init__(name, seq)

    def complement(self):  # возвращает комплементарную последовательность
        comp_seq = ''
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        for i in self.seq:
            comp_seq += comp[i]
        return comp_seq

    def reverse(self):  # разворачивает последовательность
        rev_seq = ''
        for i in self.seq:
            rev_seq = self.seq[::-1]
        return rev_seq


class Primer(DNA):

    def __init__(self, name, seq):
        super().__init__(name, seq)

    def gc_cont(self):  # считает ГЦ состав праймера и выводит процентное соотношение
        total = len(self.seq)
        c = self.seq.count("C")
        g = self.seq.count("G")
        gc_total = g + c
        gc_content = gc_total / total
        return gc_content

    def temp(self):  # считает температуру отжига праймера
        total = len(self.seq)
        a = self.seq.count("A")
        t = self.seq.count("T")
        c = self.seq.count("C")
        g = self.seq.count("G")
        temper = 22 + 1.46 * (2*(g+c) + (a+t))
        return temper


# функция для создания двух цепей ДНК
def double_dna(name, seq):
    matrix_seq = DNA(name + " matr", str(seq))  # создает матричную последовательность
    comp_seq = DNA(name, matrix_seq.complement())  # создает комплементарную последовательность
    rev_seq = DNA(name + " comp", comp_seq.reverse())  # создает обратную последовательность
    return matrix_seq, rev_seq


# обработка данных, введенных пользователем
def input_file_process(table_path, list_path):
    seqs = {}  # создали пустой словарь, чтобы добавить туда все последовательности
    with open(table_path, 'r') as table:  # открываем таблицу на чтение
        table = csv.reader(table, delimiter=';')
        next(table, None)   # пропускаем 1 строчку таблицы (названия столбцов)
        for row in table:  # пробегаемся по всей таблице, что ввел пользователь
            if row[0] not in seqs:  # проверяем наличие группы в словаре, надо вытащить из таблицы
                seqs[row[0]] = []  # если нет, то создем ключ с именем группы и к нему пустой список
            # добавляем последовательность в словарь
            with open(list_path, 'r') as list:  # открываем список на чтение
                list = SeqIO.parse(list, 'fasta')  # преобразуем fasta-файл для вычленения интересующих пар-в
                for feature in list:  # для каждой посл-ти из списка ищем совпадение имени в списке и табл
                    if feature.name in row:
                        matrix_seq, rev_seq = double_dna(str(feature.name), str(feature.seq))
                        seqs[row[0]].append(matrix_seq)
                        seqs[row[0]].append(rev_seq)

    # код для того, чтобы посмотреть словарь исходных последовательностей
    for group in seqs:
        for s in seqs[group]:
            print(group, s.name, s.seq, len(s)) # выводит группу, имя, последовательность и ее длину

    return seqs


def pre_primer_search(seq_list):
    pre_primer_set, pre_primer = set(), ''  # пустое множество для пре-праймеров и нустой препраймер
    for i in range(len(seq_list[0])):  # проходимся по всем буквам в последовательности
        letter_check = True  # изначально считаем, у всех последовательностей одинаковый i-ый нуклеотид
        for j in range(0, len(seq_list)-2, 2):  # пробегаемся по всем последовательностям
            if seq_list[j][i] != seq_list[j+2][i]:  # проверка на равенство нуклеотидов
                letter_check = False  # выявляем разницу в нуклеотидах
                if len(pre_primer) >= 16:  # если длина пре-праймера >= 16
                    pre_primer_set.add(pre_primer)  # добавляем препраймер в список
                pre_primer = ''  # обнуляем пре-праймер
                break  # обрываем цикл
        if letter_check:  # если у всех последовательностей i-ый нуклеотид одинаковый
            pre_primer += seq_list[0][i]  # записываем нуклеотид в пре-праймер
    return pre_primer_set  # возвращаем множество препраймеров


# создание forward-reverse праймеров
def forw_rev_primers(pre_primer_list):
    primers = []
    for i in pre_primer_list:  # для каждого праймера из списка
        forward = Primer('forward', i)     # создаем прямой праймер
        revers = Primer('reversed', forward.complement())  # создаем комплементарный праймер
        revers = Primer('reversed', revers.reverse())        # создаем обратный праймер
        primers.append(forward)       # добавлем в список прямые праймеры
        primers.append(revers)      # добавляем в список обратные праймеры

    return primers     # возвращаем словарь праймеров без повторов и подстрок


# фильтрация праймеров по GC составу
def gc_primer(primer_list):
    primers = []
    for i in primer_list:   # для каждого праймера из списка
        if 0.5 < i.gc_cont() < 0.6:     # если GC состав меньше 50% и больше 60%
            primers.append(i)

    # вывод праймеров
    for primer in primers:    # для каждой группы
        print(primer.name, primer.seq)       # выводит группу и отфильтрованные праймеры

    return primers   # возвращаем отфильтрованный по GC составу словарь праймеров


def primer_dict(seq_dict):
    primers_dict = {}  # создали пустой словарь для препраймеров
    for group in seq_dict:  # для каждой группы в словаре последовательностей
        pre_primers_list = list(pre_primer_search(seq_dict[group]))  # словарь группа : список пре-праймеров
        pre_primers_list.sort(key=len)  # сортируем список праймеров по возрастанию длины
        primer_list = forw_rev_primers(pre_primers_list)
        primers_dict[group] = gc_primer(primer_list)
    return primers_dict


# функция по определению, оттожется ли праймер на последовательности
def mismatch_counter(seq, primer):
    primer_annealing = False  # изначально отжиг праймера = F
    for i in range(len(seq)-len(primer)):  # чтобы "отжечь" праймер на каждом возможном участке seq
        mismatch = 0  # для подсчета побуквенных несовпадений
        for j in range(len(primer)):  # побуквенно проходимся по primer
            if primer[j] != seq[i+j]:  # побуквенно сравниваем primer и участок seq
                mismatch += 1  # считаем несовпадения
        if mismatch <= 3:  # если несовпадений было <= 3, то primer оттожется
            primer_annealing = True
            break  # обрываем цикл, так как такой праймер нам заведомо не подходит
    return primer_annealing  # возвращает T или F в зависимости от того, оттожется ли праймер на seq


# функция по созданию словаря группоспецифичных праймеров
def mismatch_sorter(seq_dict, primer_dict):
    for group in primer_dict:   # для каждой группы
        specific_primers = []  # список для специфичных праймеров
        for primer in primer_dict[group]:   # для каждого праймера из списка
            primer_annealing = False  # изначально отжиг праймера = F
            for gr in seq_dict:  # для каждой группы
                if gr != group:  # нас интересуют последовательности только других групп
                    for seq in seq_dict[gr]:  # для каждой последовательности группы
                        # проверка на отжиг праймера на последовательности
                        if mismatch_counter(seq.seq, primer.seq):
                            primer_annealing = True  # праймер отжегся
                            break  # обрываем проверку на отжиг внутри группы
                if primer_annealing:  # проверка на отжиг праймера в группе
                    break   # обрываем проверку на отжиг - этот праймер не специфичен
            if not primer_annealing:  # проверка на отжиг праймера среди всех групп
                specific_primers.append(primer)  # добавляем праймер в список специфичных праймеров
        primer_dict[group] = specific_primers  # сохраняем в словаре только группоспецифичные праймеры

    for group in primer_dict:   # для каждой группы
        for primer in primer_dict[group]:   # для каждого праймера в списке группы
            print(group, primer.name, primer.seq)

    return primer_dict  # возвращаем словарь группоспецифичных праймеров


# функция подбора пар праймеров по температуре отжига
def temp_primer_pairs(primer_list):
    pairs = []   # список пар праймеров
    for i in range(len(primer_list)):   # для каждого парймера из списка группы
        primer_1 = primer_list[i]
        for j in range(i + 1, len(primer_list)):   # проверяем каждый следующий праймер этой группы
            primer_2 = primer_list[j]
            if primer_1.name == 'forward' and primer_2.name == 'reversed':   # если могут образовать пару прмой-обратный
                if abs(primer_1.temp() - primer_2.temp()) <= 5:  # разность температуры отжига не превышает 5 по модулю
                    pairs.append(primer_1)   # добавляем прямой праймер в список
                    pairs.append(primer_2)   # добавляем обратный праймер в список
            else:
                continue   # продолжаем идти по списку, если праймеры не образуют пару прямой-обратный

    return pairs


# функция проверки взаиморасположения пары праймеров
def coord_match(seq, primer_list):
    match_list = []   # пустой список пар праймеров
    for i in range(0, len(primer_list)-1, 2):   # для каждой пары праймеров из списка группы
        primer_1, primer_2 = primer_list[i], primer_list[i+1]
        matrix_strand = seq.seq   # в переменную помещаем ее посл-ть
        comp_strand = seq.complement()   # в др переменную помещаем ей компл-ую
        primer_2_rev = Primer('reversed', primer_2.reverse())  # переворачиваем обратный
        coord1 = matrix_strand.index(primer_1.seq)  # координата forward праймера
        coord2 = comp_strand.index(primer_2_rev.seq)  # координата reversed праймера
        if coord1 != coord2:  # проверяем, что они не совпадают
            match_list.append(primer_1)  # добавляем пару праймеров
            match_list.append(primer_2)  # добавляем пару праймеров
    return match_list


# функция подбора пар праймеров
def primer_pairs(primers_dict, seq_dict):     # сортируем по температуре
    for group in primers_dict:  # для каждой группы
        temp_pairs = temp_primer_pairs(primers_dict[group])   # сохраняем пары праймеров по температуре
        primers_dict[group] = coord_match(seq_dict[group][0], temp_pairs)


        for primer in primers_dict[group]:  # для каждого праймера в списке группы
            print(group, primer.name, primer.seq, primer.temp())  # выводим: № группы, тип праймера, температуру отжига

    return primers_dict  # возвращаем словарь пар праймеров


# запись результата подбора праймеров в файл
def output_file_process(primer_pairs):
    total = []   # список списков
    groups = []   # список для записи групп
    primers_name = []   # список для записи типа праймера
    primers_seqs = []   # список для записи последовательностей праймеров
    for group in primer_pairs:   # для каждой группы из словаря пар праймеров
        for primer in primer_pairs[group]:   # для каждого праймера группы
            groups.append(group)   # записываем номер группы
            primers_name.append(primer.name)   # записываем тип праймера
            primers_seqs.append(primer.seq)   # заисываем последовательность праймера
    total.append(groups)   # добавляем списки в общий список
    total.append(primers_name)
    total.append(primers_seqs)
    export_data = zip_longest(*total, fillvalue='')
    with open('final_table.csv', 'w', encoding="ISO-8859-1", newline='') as final_table:   # создаем файл для записи
        wr = csv.writer(final_table)
        wr.writerow(("Group", "Primers", "Primer's seqs"))   # названия столбцов
        wr.writerows(export_data)   # записываем данные в файл
    return export_data


if __name__ == '__main__':
    # Пользователь вводит полные пути к файлам
    #path_to_table, path_to_list = input('Введите путь до таблицы '), input('Введите путь до списка ')
    # path_to_table, path_to_list = input('Введите путь до таблицы '), input('Введите путь до списка ')

    path_to_table = '/Users/akhvorov/Desktop/home_task/BioProject/TestPro.csv'
    path_to_list = '/Users/akhvorov/Desktop/home_task/BioProject/TestPro.fasta'
    # path_to_table = '/Users/akhvorov/Desktop/home_task/BioProject/Project.csv'
    # path_to_list = '/Users/akhvorov/Desktop/home_task/BioProject/Project_aligned.fasta'

    # path_to_table = './TestPro.csv'
    # path_to_list = './TestPro.fasta'

    # path_to_table = '/Users/dnayd/Desktop/Project/TestPro.csv'
    # path_to_list = '/Users/dnayd/Desktop/Project/TestPro.fasta'

    seq_dict = input_file_process(path_to_table, path_to_list)   # словарь исходных последовательностей
    primer_dict = primer_dict(seq_dict)
    specific_primers = mismatch_sorter(seq_dict, primer_dict)  # словарь специфичных праймеров
    matching_primers = primer_pairs(specific_primers, seq_dict)
    output = output_file_process(matching_primers)               # вывод результатов в .csv файл

