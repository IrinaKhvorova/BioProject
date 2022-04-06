import csv
from Bio import SeqIO


class Sequence:

    alphabet = []

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def __getitem__(self, item):  # возвращает по запросу символ в последовательности
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


class Primer:
    alphabet = ['A', 'T', 'G', 'C']

    def __init__(self, seq):
        self.seq = seq

    def __getitem__(self, item):  # возвращает по запросу символ в последовательности
        return self.seq[item]

    def __len__(self):  # возвращает длину последовательности
        return len(self.seq)

    def gc_cont(self):
        total = len(self.seq)
        c = self.seq.count("C")
        g = self.seq.count("G")
        gc_total = g + c
        gc_content = gc_total / float(total)
        return gc_content

    def temp(self):
        total = len(self.seq)
        a = self.seq.count("A")
        t = self.seq.count("T")
        c = self.seq.count("C")
        g = self.seq.count("G")
        temper = 22 + 1.46 * ([2*(g+c)] + (a+t))
        return temper


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
                        seqs[row[0]].append(DNA(str(feature.name), str(feature.seq)))

    # код для того, чтобы посмотреть словарь исходных последовательностей
    for key in seqs:
       for s in seqs[key]:
         print(key, s.name, s.seq, len(s)) # выводит группу, имя, последовательность и ее длину

    return seqs


# функция для разделения последовательности на k-меры
def seq_kmers(k, sequence):  # принимает длину к-мера и последовательность
    kmers = []  # создаем пустой список k-меров
    for s in range(0, len(sequence) - k + 1):  # пробегаемся по всей последовательности
        kmers.append(sequence[s:s+k])  # делаем срезы длины k и добавляем их в список
    return kmers   # возвращает списком все к-меры последовательности


# создание словаря k-меров
# в группах для каждой возможной длины праймера хранится список множеств
# в множествах находятся подстроки длины k для каждой последовательности группы
def kmers_dict(seq_dict):
    kmers_dict = {}  # создали пустой словарь k-меров
    for group in seq_dict:  # для каждой группы в словаре
        kmers_dict[group] = {}  # словарь состоит из групп
        for s in seq_dict[group]:  # для каждой последовательности в группе
            for k in range(16, 31):  # для каждой длины праймерв (задаем границы)
                if k not in kmers_dict[group]:  # проверка на наличие ключа - длины k
                    kmers_dict[group][k] = []  # если нет, то создаем ключ и пустой список для множеств
                unique_kmers = set(seq_kmers(k, s.seq))  # находим множество уникальных к-меров (длины k)
                kmers_dict[group][k].append(unique_kmers)  # добавляем в список множество с подстроками длины k
    return kmers_dict  # возвращаем словарь k-меров


# создание словаря праймеров
# оставляем из словаря k-меров только те праймеры, которые общие для всех последовательностей группы
def primer_dict(seq_dict):
    primer_dict = {}  # создаем пустой словарь праймеров
    dict = kmers_dict(seq_dict)  # получаем словарь k-меров
    for group in dict:  # для каждой группы
        primer_dict[group] = []  # словарь праймеров состоит из групп
        for k in dict[group]:  # для каждой длины праймерв в группе
            for i in range(len(dict[group][k])-1):  # пробегаемся по всем множествам k-меров в группе
                if i == 0:  # проверка на то, что это первое множество
                    k_intersection = dict[group][k][i]  # в пересечении будет первое мноество
                k_intersection.intersection_update(dict[group][k][i+1])  # находим пересечение всех множеств списка
            for primer in k_intersection:  # для каждого праймера в множестве пересейчений
                primer_dict[group].append(primer)  # добавляем праймер в словарь (группа : список праймеров)

    # код для того, чтобы посмотреть словарь праймеров
    for key in primer_dict:  # для каждой группы
        print(key, primer_dict[key]) # выводит группу и список праймеров

    return primer_dict   # возвращаем словарь праймеров


# удаляем из значений словаря праймеры, которые являются подстроками более длинных праймеров
# сортируем в порядке возрастания длины праймеров
# обновляем значения в словаре
def unique_primer(primer_dict):
    substrings_set = set()   # создаем пустое множество для добавления подстрок
    for key in primer_dict:   # для каждой группы
        primer_dict[key].sort(key=len)  # сортируем список праймеров по возрастанию длины
        for i in range(len(primer_dict[key])):   # для каждого праймера из списка
            for j in range(i + 1, len(primer_dict[key])):   # каждый следующий праймер из списка
                if primer_dict[key][i] in primer_dict[key][j]:   # проверяем входит ли в него праймер i
                    substrings_set.add(primer_dict[key][i])   # если входит, добавляем праймер i в множ-во подстрок
        substrings_list = list(substrings_set)   # полученное мн-во трансформируем в список
        # оставляем в списке только те строки(праймеры), кот-х нет в списке подстрок
        unique_primer_list = [item for item in primer_dict[key] if item not in substrings_list]
        primer_dict[key] = unique_primer_list   # перезаписываем для каждого ключа обновленный список

    # код для того, чтобы посмотреть словарь отсортированных праймеров
    for key in primer_dict:  # для каждой группы
        print(key, primer_dict[key])   # выводит группу и оставшиеся праймеры, отсорт. по длине

    return primer_dict   # возвращаем словарь праймеров без повторов и подстрок


if __name__ == '__main__':
    # Пользователь вводит полные пути к файлам
    path_to_table, path_to_list = input('Введите путь до таблицы '), input('Введите путь до списка ')

    # path_to_table = '/Users/akhvorov/Desktop/home_task/BioProject/TestPro.csv'
    # path_to_list = '/Users/akhvorov/Desktop/home_task/BioProject/TestPro.fasta'

    # Tanya's paths: table =  ./TestPro.csv  list = ./TestPro.fasta

    seq_dict = input_file_process(path_to_table, path_to_list)
    primer_dict(seq_dict)

    primers_dict = primer_dict(seq_dict)
    unique_primer(primers_dict)