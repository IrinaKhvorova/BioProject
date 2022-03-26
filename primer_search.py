# Обработка данных, введенных пользователем
def input_file_process(list_path, table_path):
    seqs = {}  # создали пустой словарь, чтобы добавить туда все последовательности
    with open(table_path) as table:  # открываем файл на чтение
        for line in table:  # пробегаемся по всей таблице, что ввел пользователь
            if '# имя группы #' not in seqs:  # проверяем наличие группы в словаре, надо вытащить из таблицы
                seqs['# имя группы #'] = []  # если нет, то создем ключ с именем группы
            # добавляем последовательность в словарь, надо вытащить из списка
            with open(list_path) as list:
                seqs['# имя группы #'].append('# сама последовательность из списка от пользователя #')


if __name__ == '__main__':
    # Здесь надо подумать, что будет вводить пользователь, думаю, что это должны быть полные пути к файлам
    path_to_table, path_to_list = input('Введите путь до таблицы '), input('Введите путь до списка ')
    input_file_process(path_to_table, path_to_list)
