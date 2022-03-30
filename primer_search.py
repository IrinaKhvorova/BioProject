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


# Обработка данных, введенных пользователем
def input_file_process(table_path, list_path):
    seqs = {}  # создали пустой словарь, чтобы добавить туда все последовательности
    with open(table_path, 'r') as table:  # открываем файл на чтение
        table = csv.reader(table, delimiter=';')
        next(table, None)   # пропускаем 1 строчку таблицы (названия столбцов)
        for row in table:  # пробегаемся по всей таблице, что ввел пользователь
            if row[0] not in seqs:  # проверяем наличие группы в словаре, надо вытащить из таблицы
                seqs[row[0]] = []  # если нет, то создем ключ с именем группы и к нему пустой список
            # добавляем последовательность в словарь
            with open(list_path, 'r') as list:
                list = SeqIO.parse(list, 'fasta')  # преобразуем fasta-файл для вычленения интересующих пар-в
                for feature in list:  # для каждой посл-ти из списка ищем совпадение имени в списке и табл
                    if feature.name in row:
                        seqs[row[0]].append(DNA(str(feature.name), str(feature.seq)))


# код для того, чтобы посмотреть наш словарь, потом уберем
    for key in seqs:
        for s in seqs[key]:
            print(key, s.name, s.seq, len(s)) # выводит группу, имя, последовательность и ее длину
    return seqs

# функция для разделения последовательности на k-меры
def seq_kmers(k, sequence):  # принимает длину к-мера и последовательность
    kmers = []
    for s in range(0, len(sequence) - k + 1):
        kmer = sequence[s:s+k]
        kmers.append(kmer)
    return kmers   # возвращает списком все к-меры последовательности

# создание словаря, где каждой возможной длине праймера в группах соответствуют пары
# (последовательность : множество подстрок длины k)
def kmers_dict(dict):
    kmers_dict = {}  # создали пустой словарь k-меров
    for group in dict:  # для каждой группы в словаре
        kmers_dict[group] = {}  # словарь состоит из групп
        for s in dict[group]:  # для каждой последовательности в группе
            for k in range(16, 31):  # для каждой длины праймерв (задаем границы)
                kmers_dict[group][k] = {}  # в группах возможные варианты длины k-меров
                unique_kmers = set(seq_kmers(k, s.seq))  # находим множество к-меров (длины k)
                kmers_dict[group][k][s] = unique_kmers  # в k хранятся (последовательность : подстроки длины k)
    print(kmers_dict)

if __name__ == '__main__':
    # Пользователь вводит полные пути к файлам
    path_to_table, path_to_list = input('Введите путь до таблицы '), input('Введите путь до списка ')
    # /Users/akhvorov/Desktop/home_task/BioProject/TestPro.csv
    # /Users/akhvorov/Desktop/home_task/BioProject/TestPro.fasta
    seq_dict = input_file_process(path_to_table, path_to_list)
    kmers_dict(seq_dict)

