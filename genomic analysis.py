from abc import ABC, abstractmethod
import re
from collections import Counter
import argparse
import sys
import os

class BIOFormat(ABC): #используем ABC - abstract class, a helper class for further inheritance
    def __init__(self, file_path, encoding=None): #пусть и кодировка вообще, VS ее определяет сама
        super().__init__()
        buf = [] #буффер, работали по аналогии
        with open(file_path, 'r', encoding=encoding) as f: #стандартное открытие файла для чтения, поэтому тут 'r'
            buf = f.readlines()
        self.data = [el.strip() for el in buf] # el = element
        self.MAX_QUALITY = 126 #значение символа '~' в ASCI
        self.MIN_QUALITY = 33 #значение символа '!' в ASCI
      
    @abstractmethod #используем декоратор, абстрактный метод
    def get_mean_gc(self):
        raise NotImplementedError('Расчет среднего GC содержимого не определен для данного типа файла') #подняли ошибку на случай
  
    @abstractmethod #используем декоратор, абстрактный метод
    def get_mean_quality_reading(self):
        raise NotImplementedError(f'Расчет среднего качества прочтений не определен для данного типа файла') #подняли ошибку на случай
  
    @abstractmethod #используем декоратор, абстрактный метод
    def get_reading_count(self):
        raise NotImplementedError(f'Расчет количества записей в файле не определен для данного типа файла') #подняли ошибку на случай
  
    @abstractmethod #используем декоратор, абстрактный метод
    def get_alignments_count(self, chromosome, start, end):
        raise NotImplementedError(f'Расчет числа выравниваний не определен для данного типа файла') #подняли ошибку на случай

    @abstractmethod #используем декоратор, абстрактный метод
    def get_filter_alignments(self, chromosome, start, end):
        raise NotImplementedError(f'Выборка выравниваний не определена для данного типа файла') #подняли ошибку на случай


class FastaFile(BIOFormat):  # наследник BIOFormat
    def __init__(self, file_path, encoding=None):  # конструктор класса
        super().__init__(file_path, encoding)  # вызов методов родительского класса
        self.sequences = list(self.split_data_sample(
            self.data))  # свойство экземпляра класса, которое хранит список всех последовательностей

    @classmethod
    def split_data_sample(cls,
                          data):  # метод, который отвечает за разбиение данных на заголовки и соответствующие им последовательности
        current_sequence = []  # пустой список создали
        current_header = None  # пустая переменная для хранения заголовка
        for line in data:  # проходим данные построчно
            if line.startswith(
                    '>'):  # если последовательность начинается с > (а если это первая строка, то ДА, начинается)
                if current_header and current_sequence:  # т.е. если они не None
                    yield [current_header, ''.join(
                        current_sequence)]  # возвращает текущую пару (заголовок и последовательность) в виде списка, (yield возвращает значения по-одному и запоминает где остановилась, т.к. это функция-генератор)
                current_header = line  # обновляем текущий заголовок
                current_sequence = []  # обнуляем последовательность
            else:
                current_sequence.append(
                    line.strip())  # если строка - не заголовок: добавляем ее в текущую последовательность, убираем пробелы функцией strip()

        if current_header and current_sequence:  # проверка для возврата последней последовательности, если она существует
            yield [current_header, ''.join(current_sequence)]

    def get_mean_gc(self):  # метод для расчета среднего GC-контента
        data_count = {}
        for header, sequence in self.sequences:
            gc_content = (sequence.count('G') + sequence.count('C')) / len(
                sequence)  # считаем G,C в последовательности, затем делим на длину последовательности
            data_count[
                header] = gc_content  # в словаре data_count заголовок последовательности (header) используется как ключ, а рассчитанное содержание GC (gc_content) — как значение

        return [f'{key}: {value}' for key, value in
                data_count.items()]  # цикл, который проходит по всем элементам словаря, возвращает 'заголовок: последовательность'

    def get_mean_quality_reading(self):  # для моего класса не нужно
        return 0

    def get_reading_count(self):  # возвращает длину последовательности
        return len(self.sequences)

    def get_alignments_count(self, chromosome=None, start=None, end=None):  # для моего класса не нужно
        return 0

    def get_filter_alignments(self, chromosome=None, start=None, end=None):  # для моего класса не нужно
        return []

class FastqFile(FastaFile): #потомок Фасты. Наследование классов
    def __init__(self, file_path, encoding=None): #все то же, что иу Фасты. Только для удобство данные будут лежать в data_fastq в виде списка по 4 строки каждый абзац
        super().__init__(file_path, encoding) #encoding = это кодировка файла, она разная у каждого тут
        self.data_fastq = list(self.split_data_sample(self.data)) #делим методом, в отдельнную переменную пихаю

    @classmethod
    def split_data_sample(self, data): #метод класса юзаю, поэтому тут и декоратор 
        for i in range(0, len(data), 4): #4 строки, так как файл весь из абзацев по 4 строки нуждо делить
            yield data[i:i+4] #генератор, возвращает значения один раз и не хранит в памяти, удобно.

    def get_mean_gc(self): # гуанин-цитозин контент подсчитывает и усредняет
        data_count = {} #в словарь запихну
        for i, value in enumerate(self.data_fastq): #enumerate() функция вообще возвращает кортежи по типу (а, б), i=мой "номер", а value="содержимое"
            data_count[i] = (value[1].count('G') + value[1].count('C'))/len(value[1]) #для value мне НУЖЕН ИНДЕКС 1 !!!
        return [f'{i}: {value}' for i, value in data_count.items()] #возвращает контент для каждой строки для красоты, но выглядит все равно огромной махиной
  
    def get_mean_quality_reading(self): #среднее качество прочтений, map() - даст возможность применить функцию к переменным (итерируемой)
         # порядковый номер: качество. Пробегает по элементу каждому в data_fastq, подсчитывая качество. Качество прочтения лежит в 4 строке каждого "абзаца"
        data_quality = {i+1: (sum(ord(q) - self.MIN_QUALITY for q in value[3]) / len(value[3])) / self.MAX_QUALITY for i, value in enumerate(self.data_fastq)} 
        #нормировка качествоа секвенироования, от 0 до 100
        #тут решила идти через индексы, можно же 4 строку через value[3], начинается же с 0
        return data_quality
  
    def get_reading_count(self):
        return len(self.data_fastq) #решила добавить еще вот это, кол-во прочтений

    def get_alignments_count(self, chromosome=None, start=None, end=None): #для моего класса это как бы не нужно, добавила параметры, по умолчанию пусто
        return 0

    def get_filter_alignments(self, chromosome=None, start=None, end=None): #для моего класса это не нужно снова, добавила параметры, по умолчанию пусто
        return []
  
class VcfFile(BIOFormat):#класс для работы с VCF файлами, наследует функциональность от BIOFormat
    def __init__(self, file_path, encoding=None):#file_path-путь к VCF файлу; encoding - кодировка файла, по умолчанию None.
        super().__init__(file_path, encoding)#инициализация родительского класса
        vsf_data = list(filter(lambda x: x[0] != '#', self.data))#фильтруем данные, исключая строки, начинающиеся с '#'тк это комментарии
        self.vsf_data = [row.split('\t') for row in vsf_data]#разбиваем оставшиеся строки по символу табуляции и сохраняем в атрибуте vsf_data

    def get_mean_gc(self):#среднее содержание GC в VCF
        gc_row = [row for row in self.vsf_data if row[3] + row[4] == 'GC']#фильтруем строки, где 4-й и 5-й элементы 'GC'
        return len(gc_row) / len(self.vsf_data)#отношение количества GC строк к общему количеству строк

    def get_alignments_count(self, chromosome, start, end):#количество выравниваний на заданном хромосоме в указанном диапазоне
        alignments = [
            row for row in self.vsf_data 
            if row[0] == chromosome and int(row[1]) >= start and int(row[1]) <= end#фильтруем данные по заданным хромосоме и диапазону
        ]
        return len(alignments)#количество найденных выравниваний

    def get_mean_quality_reading(self):#среднее качество чтения.
        qualities = [float(row[5]) for row in self.vsf_data if row[5] != '.']#row[5] -предполагаемый индекс, где инфa о качестве чтения, Добавляем в список значения, не равные '.'
        return sum(qualities) / len(qualities) if qualities else 0#если значения корректны-возвращаем среднее значение; если нет-возвращаем 0

    def get_filter_alignments(self, chromosome, start, end):#фильтрует выравнивания по заданным критериям и возвращает их в виде строк
        alignments = [
            row for row in self.vsf_data 
            if row[0] == chromosome and int(row[1]) >= start and int(row[1]) <= end
        ]#фильтруем данные по заданным хромосоме и диапазону
        return [' '.join(row) for row in alignments]#список строк, элементы каждой строки соединены пробелом

    def get_reading_count(self):#количество чтений для каждого хромосомы
        chromosomes = Counter([el[0] for el in self.vsf_data])#количество вхождений каждой хромосомы
        return [f'{key}: {value}' for key, value in chromosomes.items()]#форматируем результат для удобного отображения как словарь
  
class SamFile(BIOFormat):# наследник BIOFormat
    def __init__(self, file_path, encoding=None):# конструктор класса
        super().__init__(file_path, encoding)# вызов методов родительского класса
        sam_data = list(filter(lambda x: x[0]!='@', self.data))[1:]#lambda-фильтрует строки, чтобы исключить заголовки 
        self.sam_data = [row.split('\t') for row in sam_data]#содержит данные SAM, разбивая каждую строку на элементы 

    def get_mean_gc(self):# метод для расчета среднего GC-контента
        data_count = {}
        for i, value in enumerate(self.sam_data):#enumerate() функция возвращает кортежи
            data_count[i] = (value[9].count('G') + value[9].count('C'))/len(value[9])# Подсчет количества 'G' и 'C' в 10-м столбце
        return [f'{key}: {value}' for key, value in data_count.items()]# цикл, который проходит по всем элементам словаря, возвращает 'заголовок: последовательность'
        # в дата каунт cчитаем количество G и C, делим на длину последовательности
    
    def get_alignments_count(self, chromosome, start, end):# Подсчет количества выравниваний на заданном хромосоме в указанном диапазоне
        alignments = [#Создам список 
            row for row in self.sam_data #Это конструкция "генератора списка- oна перебирает все строки
            if row[2] == chromosome and int(row[3]) >= start and int(row[3]) <= end
        ]
        return len(alignments)
    #row проверяем равен ли третий элемент строки (индекс 2) значению переменной (что мы ищем выравнивания, относящиеся к конкретной хромосоме)
    #init-start проверяеь находится ли позиция выравнивания (четвертый элемент строки, индекс 3) в данном диапазоне 
    #int(row[3]) преобразует значение позиции из строки в целое число
  
    def get_mean_quality_reading(self):#среднее качество прочтений
         # порядковый номер: качество. Пробегает по элементу каждому в data_fastq, подсчитывая качество. Качество прочтения лежит в 4 строке каждого "абзаца"
        data_quality = {#Создам словарь data_quality, в котором будут храниться средние значения качества прочтений
            i+1: (sum(ord(q) - self.MIN_QUALITY for q in el[10]) / len(el[10])) / self.MAX_QUALITY #Ключ словаря — это порядковый номер рида
            for i, el in enumerate(self.sam_data)#перебрает каждую строку
        }
        return data_quality
    #Для каждого символа качества q в 11-м столбце (индекс 10) рида el вычисляется его числовое значение с помощью функции ord-вычитем для норм вида-делим на колво симв в 11 столбце=ср знач
    #делим self для норм вида в диапазон 0-1
  
    def get_filter_alignments(self, chromosome, start, end):#возвращает отфильтрованные выравнивания в виде строк
    #строка соотв данной хромосоме и диапазону
        alignments = [
            row for row in self.sam_data 
            if row[2] == chromosome and int(row[3]) >= start and int(row[3]) <= end
        ]
        return [' '.join(row) for row in alignments]
    #Метод join() объединяет все элементы списка row в одну строку, вставляя между ними пробелы.
  
    def get_reading_count(self):# возвращает длину последовательности
        chromosomes = Counter([el[2] for el in self.sam_data])# создает объект Counter из модуля collections, который подсчитывает количество вхождений каждого элемента
        return [f'{key}: {value}' for key, value in chromosomes.items()]#chromosomes.items() возвращает пары (ключ, значение) для кажд уник хромосомы и еe количества.
    #f'{key}: {value}' создает строку, где key — это имя хромосомы, а value — количество её вхождений.
    


def createParser(): #некий navigation
    #к примеру, тут есть навигация по существующим командам !!!
    # 'GC_content', 'Quality', 'SliceFile', 'AlignmentsCount', 'SliceCount'
    parser = argparse.ArgumentParser(
        prog = 'main.py', # Указываем имя программы
        description = 'Программа для работы с биологическими файлами форматов fasta, fastq, sam, vsf' #описание программы
    )
    # + обяз аргумент 'function'(какую функцию выполнить)
    parser.add_argument(
        'function', #имя аргумента
        help="Функция, применяемая к файлу", #описание аргумента,что видим в справке
        choices=['GC_content', 'Quality', 'SliceFile', 'AlignmentsCount', 'SliceCount'], 
        metavar = 'function_name' #какое имя аргумента видим в выводе справки
    )
    #+ обяз аргумент 'file_path'(какой путь к входному файлу)
    parser.add_argument('file_path', help="Путь к входному файлу", metavar = 'input_file_path')#описание и имя аргумента
    # + необяз аргумент '-o' или '--output'(путь к выходному файлу)
    parser.add_argument('-o','--output', help="Путь к выходному файлу", metavar = 'output_file_path')
    # + необяз аргумент '-r' или '--range'(какой диапазон выравнивания)
    parser.add_argument('-r', '--range', help="Диапазон выравнивания", metavar='aligments_range')
    return parser

def get_example_files_path(): #функция специально для получения наших путей к примерам-файлам. Долждны быть подгружены !!!
    base_path = os.path.join(os.path.dirname(__file__), 'example_files') #определяем базовый путь к папке с примерами файлов
    
    # Возвращаем словарь с путями к примерам файлов
    return {
        'fasta': os.path.join(base_path, 'example.fasta'),
        'fastq': os.path.join(base_path, 'example.fastq'),
        'sam': os.path.join(base_path, 'example.sam'),
        'vcf': os.path.join(base_path, 'example.vcf')
    }

if __name__ == '__main__': #стандартная наша любимая фраза
    try:
        parser = createParser()#создаем парсер командной строки для обработки аргументов
        namespace = parser.parse_args(sys.argv[1:]) #нет никакого параметра. передадим в него все параметры из sys.argv кроме нулевого(содержит имя нашей программы)
        #sys.argv содержит список параметров, переданных программе через командную строку

        #определяем тип файла на основе его расширения и создаем соответствующий объект файла
        if namespace.file_path.endswith('.fasta'): #определяет, заканчивается ли данная строка определённым суффиксом- расширением фасты. далее то же самое для других файлов
            file = FastaFile(namespace.file_path) #создает объект класса `FastaFile`, передавая в конструктор путь к файлу
        elif namespace.file_path.endswith('.fastq'):
            file = FastqFile(namespace.file_path, encoding='utf-16LE') #для fastq с заданной кодировкой
        elif namespace.file_path.endswith('.sam'):
            file = SamFile(namespace.file_path, encoding='utf-16LE')#для sam с заданной кодировкой
        elif namespace.file_path.endswith('.vcf'):
            file = VcfFile(namespace.file_path)#для vcf

#в зависимости от указанной функции выполняем соответствующий метод объекта файла
        if namespace.function == 'GC_content':
            res = file.get_mean_gc()#получаем среднее содержание GC

        elif namespace.function == 'Quality':
            res = file.get_mean_quality_reading()#получаем среднее значение качества чтения

        elif namespace.function == 'SliceFile':
            region = namespace.range#получаем диапазон из аргументов
            chromosome, positions = region.split(':')#разделяем хромосому и позиции
            start, end = map(int, positions.split('-'))#преобразуем позиции в целые числа
            res = file.get_filter_alignments(chromosome, start, end)#получаем выравнивания в заданном диапазоне

        elif namespace.function == 'AlignmentsCount':#получаем общее количество прочтений
            res = file.get_reading_count()

        elif namespace.function == 'SliceCount': #то же, что для 'SliceFile'
            region = namespace.range
            chromosome, positions = region.split(':')
            start, end = map(int, positions.split('-'))
            res = file.get_alignments_count(chromosome, start, end)# Получаем количество выравниваний в заданном диапазоне

        if namespace.output: #это обращение к атрибуту с именем `output` внутри объекта `namespace`
            with open(namespace.output, 'w') as f:#открываем файл для записи
                for line in res:#проходим по результатам
                    f.write(f"{line}\n")#записываем каждую строку результата в файл
        else:
            print(res)#если выходной файл не указан, выводим результат на экран

    except Exception as e:#обрабатываем любые исключения
        print(e)#выводим сообщение об ошибке
