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

    @property
    def len(self):
        return len(self.seq)  # возвращает длину последовательности

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
def input_file_process(list_path, table_path):
    seqs = {}  # создали пустой словарь, чтобы добавить туда все последовательности
    with open(table_path) as table:  # открываем файл на чтение
        for line in table:  # пробегаемся по всей таблице, что ввел пользователь
            if '# имя группы #' not in seqs:  # проверяем наличие группы в словаре, надо вытащить из таблицы
                seqs['# имя группы #'] = []  # если нет, то создем ключ с именем группы и к нему пустой список
            # добавляем последовательность в словарь, надо вытащить из списка
            with open(list_path) as list:
                seqs['# имя группы #'].append(DNA('# имя последовательности #', '# сама последовательность #'))


if __name__ == '__main__':
    # Здесь надо подумать, что будет вводить пользователь, думаю, что это должны быть полные пути к файлам
    path_to_table, path_to_list = input('Введите путь до таблицы '), input('Введите путь до списка ')
    input_file_process(path_to_table, path_to_list)
