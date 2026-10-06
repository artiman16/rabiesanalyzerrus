import os
from os import path
import time
import shutil
import re
import argparse, sys
import glob
import json
import subprocess

version = 1.0 # Убрано влияние стоп-кодонов на классификацию, но указывается также в валидации

try: #Проверяем наличие графической библиотеки
        import dearpygui.dearpygui as dpg
except ImportError as e:
        print(f"'{e.name}' не установлен. Вы можете использовать только CLI режим с параметрами")

try: #Проверяем наличие пакетов для анализа
        import pandas as pd
        import numpy as np
        from Bio import AlignIO
        from io import StringIO
        from tqdm import tqdm
except ImportError as e:
        print(f"Ошибка: '{e.name}' не установлен. Установите пакеты с помощью предоставленного environment.yml.")
        sys.exit(1)

isGui = False
if shutil.which("mafft") is None: #Проверяем наличие mafft для выравнивания
        print("Ошибка: mafft не найден! Работа программы невозможна", file=sys.stderr)
        sys.exit(1)

def smart_sleep(seconds=1): #Если работаем в консольном варианте, ускоряем анализ
        if isGui:
                time.sleep(seconds)
                
#Обозначаем требуемые параметры для работы консольной версии
parser=argparse.ArgumentParser(description=f"RabiesAnalyzerRus v.{version}  \n - Program for classification of Lyssavirus rabies genomes by N and G gene's SNP's (Author: A.A. Gerasimenko, Rostov-on-Don Antiplague Institute, Russia)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-f", help="Input file", dest='file')
parser.add_argument("-d", help="Input directory", dest='dir')
parser.add_argument("-o", help="Output directory", dest='out', type=str, default='Results/')
parser.add_argument("-p", help="Output format (tsv or json)", dest='format', type=str, default='tsv')

args=parser.parse_args()

os.makedirs(args.out, exist_ok=True)

path = os.path.abspath(os.path.dirname(__file__))
os.chdir(path)

# Глобальные переменные
thresholdN = 11
thresholdG = 10

total_genomes = 0
progress_step = 0.0
current_progress = 0.0

our_class = "ClustersN.csv"
our_class_G = "ClustersG.csv"
troup_class = "Clusters_Troupin.csv"

nStart = 70
nEnd = 1422
gStart = 3317
gEnd = 4891

nLen = int(((nEnd - nStart + 1)/3)-1)
gLen = int(((gEnd - gStart + 1)/3)-1)

table_row = {'Sample name': '', 
             'N gene': '', 
             'G gene' : '',
             'N classification' : '', 
             'G classification' : '', 
             'By Troupin C.' : '',
             'MutationsN': '',
             'MutationsG': '',
             }

#список корявых кодонов
amb = ['---', 'a--', 't--', 'g--', 'c--', 
       'aa-', 'at-', 'ag-', 'ac-',
       'ta-', 'tt-', 'tg-', 'tc-',
       'ga-', 'gt-', 'gg-', 'gc-',
       'ca-', 'ct-', 'cg-', 'cc-']

def printing(gui, element, text, pbar=None): #Функция вывода информации в графике или консоли
        if gui == True:
                dpg.set_value(element, text)
        else: 
                if pbar: 
                        pbar.write(text)
                else: 
                        print(text)

def comp(list1, list2): #Функция поиска корявых кодонов
        for val in list1:
                if val in list2:
                        return True
        return False

def deal_ambigous(seq): #Функция замены корявых кодонов на временные нули, убираем делеции в середине генома
        for n in range(0, len(seq)-1):
                if (comp(amb,seq[n]) and seq[n+1] == '---'): 
                        seq[n] = "000"
                        seq[n+1] = "000"
                new_str = "".join(seq)
                new_str= new_str.replace("---", "000")
                new_seq = new_str.replace("-", "")

        return new_seq

def translate_dna(sequence): #Функция трансляции нуклеотидных кодонов в аминокислоты

        codontable = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                      'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                      'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                      'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                      'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                      'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                      'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                      'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                      'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                      'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                      'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                      'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                      'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                      'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                      'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
                      'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
                      '---': '-',
                      }

        seq = sequence.upper()
        prot = []

        for n in range(0, len(seq), 3):
                if seq[n:n + 3] in codontable:
                        residue = codontable[seq[n:n + 3]]
                else:
                        residue = "X"

                prot.append(residue)

        return "".join(prot)

def check_list(lst): #Функция проверки, есть ли у нас вообще ген в сиквенсе
        if not lst:
                return True  # Обработка пустого списка (можно заменить на False при необходимости)
        first = lst[0]
        return all(sublist == first for sublist in lst)

def check_gene(ref, genome, gene, gene_start, gene_end, expected_len): #Основная функция проверка качества
        proverka = True
        #Разбиваем "присланный" ген на кодоны
        res = [genome[i:i+3] for i in range(gene_start, gene_end, 3)]

        if check_list(res):
                res = [ref[i:i+3] for i in range(gene_start, gene_end, 3)]
                proverka = False

        #Заменяем пропуски (кривые кодоны) на 000
        sequ2 = deal_ambigous(res)
        res2 = [sequ2[i:i+3] for i in range(0, len(sequ2), 3)]

        #Транслируем "исправленную" последовательность
        itog = (translate_dna(deal_ambigous(res2)))

        stroka = list(itog)

        if any(item == "*" for item in stroka):
                new_itog = itog.split("*")
                pos = len(new_itog[0])

                #Заменяем "транслированные" аминокислоты после стоп-кодона иксами
                for j in range(len(new_itog[0])+1,len(itog)):
                        stroka[j] = "X"
        else:
                pos = len(stroka)-1
                stroka[len(stroka)-1] = "."

        #Печатаем получившийся протеин и его длину
        protein = "".join(stroka)
        prot = re.split('[.*]', protein)[0]

        #Считаем, сколько у нас неотсеквенировано в начале и конце генома
        begin = 0
        end = 0

        for i in range(0, len(stroka)):
                if stroka[i] == "X":
                        begin+=1
                else:
                        break

        for k in range(len(stroka)-2,0,-1):
                if stroka[k] == "X":
                        end+=1
                else:
                        break

        if (pos != len(protein)-1):
                warn = (f"Стоп-кодон в {pos+1}/{expected_len} позиции белка?")
        else: 
                warn = ("Полный")

        if proverka == False:
                dict1 = {'Gene': gene, 'Reference protein length': nLen, 'Sample protein length': 0, 
                         'Stop codon': 'Ген не найден', 'N in the begin': '-', 'N in the end': '-', 'Protein': 'Protein not found'}
        else:
                dict1 = {'Gene': gene, 'Reference protein length': gLen, 'Sample protein length': len(prot), 
                         'Stop codon': warn, 'N in the begin': begin, 'N in the end': end, 'Protein': prot}

        return dict1

def check_species(ref, probe, gene, gene_start, gene_end): #Функция проверки, бешенство ли вообще перед нами
        ref_n = ref[gene_start:gene_end].upper()
        probe_n = probe[gene_start:gene_end].upper()	
        matches = sum(1 for a, b in zip(ref_n, probe_n) if a == b and a != '-')
        identity = matches / len(ref_n) if len(ref_n) > 0 else 0
        identity = round(identity, 2)
        if identity < 0.8: 
                val_error = f"Низкая идентичность в {gene}-гене ({identity*100}% совпадения с референсом)"
        else:
                val_error = 'Валидация пройдена'
        return val_error

def snp_distance(mut_a, mut_b):
        # Преобразуем в словари: позиция → нуклеотид
        def to_dict(mut_list):
                d = {}
                for m in mut_list:
                        # Извлекаем позицию и alt
                        i = 1
                        while i < len(m) and m[i].isdigit():
                                i += 1
                        pos = int(m[1:i])
                        alt = m[i:]
                        d[pos] = alt
                return d

        dict_a = to_dict(mut_a)
        dict_b = to_dict(mut_b)

        all_pos = set(dict_a.keys()) | set(dict_b.keys())
        dist = 0
        for p in all_pos:
                base_a = dict_a.get(p)
                base_b = dict_b.get(p)
                if base_a != base_b:  # один None, другой строка → различие; или разные строки
                        dist += 1
        return dist

def countRes(row, table, nameTable, mutColumn): #Находим ближайший похожий штамм
        # Получаем мутации пробы
        mutsProbe = row[mutColumn]  

        # Читаем справочную таблицу
        ref_table = pd.read_csv(table, sep='\t')

        # Преобразуем колонку 'mutations' в списки
        mutsRef = []
        for mut_str in ref_table[mutColumn]:
                if pd.isna(mut_str) or mut_str == '':
                        mut_list = []
                else:
                        mut_list = [m.strip().upper() for m in str(mut_str).split(',') if m.strip()]
                mutsRef.append(mut_list)

        # Считаем SNP-дистанцию ОДИН РАЗ для всех строк
        ref_table['snp_distance'] = [snp_distance(mutsProbe, ref_mut) for ref_mut in mutsRef]

        # Сортируем по дистанции
        ref_table.sort_values("snp_distance", ascending=True, inplace=True)
            
        #Возвращаем данные первой (ближайшей) строки
        grupa = str(ref_table['cluster'].iloc[0])
        
        return grupa

def analysis(row): #Определяем, к какой группе принадлежит штамм
        groupOurN = countRes(row, our_class, "Ngene_Gerasimenko", "MutationsN")
        groupOurG = countRes(row, our_class_G, "Ggene_Gerasimenko", "MutationsG")
        groupTroupin = countRes(row, troup_class, "Troupin", "mutations")
        full = [groupOurN, groupOurG, groupTroupin]
        
        return full 

def resultat(row): #Расшифровываем, к какой группе относится штамм
        resultat = analysis(row)
        smart_sleep(1)
        return resultat

def json_format(df): #Функция записи таблицы в формате JSON
        df.columns = [re.sub(r'[\t\n\r\s]+', ' ', str(col).strip()).strip() for col in df.columns]
        df = df.loc[:, df.columns.notna() & (df.columns != '')]

        records = df.to_dict(orient='records')

        output_file = 'AboutProbes.json'
        with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(records, f, indent=2, ensure_ascii=False)

        return output_file

def table_format(df): #GUI: Функция сохранения выходного файла в json формате 
        item = dpg.get_value("table_format_selector")
        print(f'Вы выбрали формат {item}')
        if item == 'json': 
                res_file = json_format(df)               
                print(f'File will be saved in json format')
        else:
                res_file = 'AboutProbes.tsv'
                df.to_csv(res_file, sep='\t', index=False)
                print(f'Таблица будет сохранена в tsv формате')
        return res_file

def make_muts(df0, start, end):
        arr = []
        if (df0.empty): reps = ''
        else:
                df0 = df0[df0.POS.values >= start]
                df0 = df0[df0.POS.values <= end]
                df0.POS = df0.POS - start
                x = df0.to_string(header=False, index=False, index_names=False).split('\n')
                reps = [''.join(ele.split()) for ele in x]
                arr.append(reps)        

        if arr and len(arr[0]) > 0:
                mutations_list = [mut.strip().upper() for mut in arr[0] if mut.strip()]
        else:
                mutations_list = []
                
        return mutations_list
                
def push(): #Основная функция анализа данных
        global file_list
        ref_file = str(os.path.abspath('rabies_ref.fasta'))

        if isGui == False:
                # Создаём кортеж состояний: (file_дан?, dir_дан?)
                file_state = (args.file is not None, args.dir is not None)

                if file_state == (True, False):      # Только файл
                        file_list = [args.file]
                        print(f"📄 Один файл: {args.file}")

                elif file_state == (False, True):    # Только папка
                        file_list = glob.glob(os.path.join(args.dir, "*.fa*"))
                        print(f"📁 Папка: {args.dir} → {len(file_list)} файлов")

                elif file_state == (True, True):     # И файл, И папка
                        full_path = os.path.join(args.dir, os.path.basename(args.file))
                        if os.path.exists(full_path):
                                file_list = [full_path]
                                print(f"📄 Файл из папки: {full_path}")
                        else:
                                print(f"❌ Файл {args.file} не найден в {args.dir}")
                                exit(1)

                else:  # (False, False) — ничего не задано
                        print("❌ Введите параметры для работы программы!")
                        exit(1)

        count_mistakes = 0
        all_rows = []
        
        with tqdm(total=len(file_list), desc="Анализ геномов", 
                  bar_format="{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}") as pbar:        

                for sam in range(len(file_list)):
                        if isGui: dpg.show_item("loader")
                        in_file = str(file_list[sam])
                        name = os.path.splitext(os.path.basename(file_list[sam]))[0]		
        
                        printing(isGui, "Analysis",f'Выравнивание {name}...', pbar=pbar)
        
                        cmd = [
                                "mafft",
                            "--quiet",
                            "--reorder",
                            "--6merpair",
                            "--keeplength",
                            "--addfragments",
                            in_file,
                            ref_file
                        ]
        
                        # Запуск
                        result = subprocess.run(
                                cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True
                        )
        
                        # Обработка ошибок
                        if result.returncode != 0:
                                print("MAFFT stderr:", result.stderr)
                                raise RuntimeError(f"MAFFT failed with code {result.returncode}")
        
                        # Преобразование в объект Biopython
                        alignment = AlignIO.read(StringIO(result.stdout), "fasta") 
        
                        smart_sleep(1)
                        printing(isGui, "Analysis",f'Выравнивание {name} завершено!', pbar=pbar)
                        smart_sleep(1)
        
                        printing(isGui, "Analysis",f'Валидация {name}...', pbar=pbar)
        
                        ref = str(alignment[0].seq)
                        sequ = str(alignment[1].seq)
        
                        current_row = table_row.copy()
        
                        val_n = check_species(ref, sequ, 'N', nStart, nEnd)
                        val_g = check_species(ref, sequ, 'G', gStart, gEnd)
        
                        if (val_n or val_g) != 'Валидация пройдена':
                                count_mistakes += 1 
                                printing(isGui, "Analysis",f'Валидация {name} не пройдена!\n{val_n}\n{val_g}', pbar=pbar)
                                smart_sleep(3)
        
                        about_N = check_gene(ref, sequ, "N gene", nStart, nEnd, nLen)
                        about_G = check_gene(ref, sequ, "G gene", gStart, gEnd, gLen)
        
                        qc_result = qcForTable(name, about_N['Stop codon'],about_G['Stop codon'])
                        smart_sleep(1)
        
                        if (val_n or val_g) != 'Валидация пройдена':
                                current_row['Sample name'] = name
                                current_row['N gene'] = 'Валидация не пройдена'
                                current_row['G gene'] = 'Валидация не пройдена'
                        
                        else:
                                current_row['Sample name'] = name
                                current_row['N gene'] = about_N['Stop codon']
                                current_row['G gene'] = about_G['Stop codon']                    
        
                        if (val_n or val_g) != 'Валидация пройдена': 
                                all_rows.append(current_row)
                                continue              
        
                        printing(isGui, "Analysis",f'Валидация {name} завершена!', pbar=pbar)
                        smart_sleep(1)
        
                        printing(isGui, "Analysis",f'Поиск SNP в {name}...', pbar=pbar)
                        smart_sleep(1)
        
                        rows = []
                        rows_list = []
                        for i in range(0, len(alignment[0])):
                                if alignment[1][i] != alignment[0][i]: 
                                        dict1 = {'REF': alignment[0][i], 'POS': i, 'ALT': alignment[1][i]} 
                                        rows_list.append(dict1)
                        df0 = pd.DataFrame(rows_list) 
                        df0.ALT = df0.ALT.replace('-', 'del')
                        df0.REF = df0.REF.replace('-', 'ins')
                        k = 0
                        ins = 0
        
                        # Проходим по строкам DataFrame
                        for ii in range(len(df0)):
                                if df0.iloc[ii, df0.columns.get_loc('ALT')] == "del":
                                        k += 1
        
                                if df0.iloc[ii, df0.columns.get_loc('REF')] == "ins":
                                        pos_col = df0.columns.get_loc('POS')
                                        df0.iloc[ii, pos_col] = df0.iloc[ii, pos_col] - k + 2
                                        ins += 1
        
                                if df0.iloc[ii, df0.columns.get_loc('ALT')] != "del":
                                        pos_col = df0.columns.get_loc('POS')
                                        df0.iloc[ii, pos_col] = df0.iloc[ii, pos_col] - ins + 1
                                               
                        # Сохраняем как список, НЕ как строку в списке!
                        current_row['mutations'] = make_muts(df0, nStart, nEnd)
                        current_row['MutationsN'] = make_muts(df0, nStart, nEnd)
                        current_row['MutationsG'] = make_muts(df0, gStart, gEnd)
                         
                        printing(isGui, "Analysis",f'Поиск SNP в {name} завершен!', pbar=pbar)             
                        smart_sleep(1)
        
                        printing(isGui, "Analysis",f'Определение группы для {name}...', pbar=pbar)
                        smart_sleep(1)
                        
                        valki = resultat(current_row)
                        
                        #n_is_bad = (val_n != 'Валидация пройдена') or (about_N['Stop codon'] != 'Полный')
                        #g_is_bad = (val_g != 'Валидация пройдена') or (about_G['Stop codon'] != 'Полный')
                
                        #current_row['N classification'] = "\u00D7" if (n_is_bad == True or n_is_bad == False) else valki[0]
                        #current_row['G classification'] = "\u00D7" if (g_is_bad == True or g_is_bad == False) else valki[1]
                        #current_row['By Troupin C.'] = "\u00D7" if (n_is_bad == True or n_is_bad == False) else valki[2]
                        
                        n_is_bad = (val_n != 'Валидация пройдена')
                        g_is_bad = (val_g != 'Валидация пройдена')
                        
                        current_row['N classification'] = "\u00D7" if n_is_bad else valki[0]
                        current_row['G classification'] = "\u00D7" if g_is_bad else valki[1]
                        current_row['By Troupin C.'] = "\u00D7" if n_is_bad else valki[2]
                        
                        printing(isGui, "Analysis",f'Группа для {name} определена!', pbar=pbar)

                        all_rows.append(current_row) # СКЛАДЫВАЕМ ВСЕ НАШИ ДАННЫЕ В СТРОКУ ДАТАФРЕЙМА

                        smart_sleep(1)

                        if isGui:
                                with dpg.table_row(parent="quality_table"):
                                        for idx, value in enumerate(qc_result):
                                                if idx == 0:	
                                                        dpg.add_text(value)	
                                                else:
                                                        color = colorByvalue(value)
                                                        dpg.add_text(value, color=color)

                                vals = [current_row['Sample name'], 
                                        current_row['N classification'], 
                                        current_row['G classification'],                                
                                        current_row['By Troupin C.']]

                                with dpg.table_row(parent="class_table"):
                                        for value in vals:
                                                dpg.add_text(value)			

                                global total_genomes, current_progress, progress_step, sample_name

                                current_progress += progress_step
                                dpg.set_value("progress_bar", current_progress)
                                dpg.set_value("progress_text", f"Анализ: {int(current_progress * 100)}%")

                                dpg.hide_item("loader")
                                
                        pbar.update(1)

        df = pd.DataFrame(all_rows)
        
        df['mutations'] = df['mutations'].apply(lambda x: ', '.join(x) if isinstance(x, list) else '')
        df['MutationsN'] = df['MutationsN'].apply(lambda x: ', '.join(x) if isinstance(x, list) else '') 
        df['MutationsG'] = df['MutationsG'].apply(lambda x: ', '.join(x) if isinstance(x, list) else '') 

        df = df.drop(columns=['mutations', 'MutationsN', 'MutationsG'])
        
        df.loc[df['N gene'] == 'Валидация не пройдена', ['N classification', 'By Troupin C.']] = "\u00D7"
        df.loc[df['G gene'] == 'Валидация не пройдена', ['G classification']] = "\u00D7"
        
        df = df.rename(columns={'Sample name': 'Штамм', 'N gene': 'N ген', 'G gene': 'G ген', 
                           'N classification': 'Классификация по N гену' , 'G classification': 'Классификация по G гену', 
                           'By Troupin C.': 'Международная классификация по N гену'})
        
        if isGui: 
                res_file = table_format(df)
                shutil.move(res_file, f'Results/{res_file}')

        if isGui == False and args.format == 'json': 
                res_file = json_format(df)
                shutil.move(res_file, f'{args.out}/{res_file}')

        if isGui == False and args.format == 'tsv':
                df.to_csv(f'{args.out}/AboutProbes.tsv', index=False)

        if (count_mistakes) == 0:	
                printing(isGui, "Analysis",'АНАЛИЗ ЗАВЕРШЕН!', pbar=pbar)
        else: printing(isGui, "Analysis",f'АНАЛИЗ ЗАВЕРШЕН! Валидацию не прошли {count_mistakes} посл.', pbar=pbar)

def qcForTable(name, aboutN, aboutG): #GUI: Заполняем таблицу с информацией о качестве генома
        qual = ('Пройдена' if aboutN == 'Полный' and aboutG == 'Полный' else
        'Ошибки в N и G гене' if aboutN != 'Полный' and aboutG != 'Полный' else
        'Ошибка в N или G гене')

        values = [name, aboutN, aboutG, qual]
        return values

def colorByvalue(value): #GUI: Раскрашиваем ячейки "качества" генов
        return ([0, 255, 0, 255] if value in ("Полный", "Пройдена") 
                else [255, 165, 0, 255] if value == 'Ошибка в N или G гене' 
                else [255, 0, 0, 255])

def run_gui(): #Функция запуска графического режима      	
        dpg.create_context()
        def set_light_theme():
                with dpg.theme() as light_theme:
                        with dpg.theme_component(dpg.mvAll):
                                # Настраиваем цвета (пример)
                                dpg.add_theme_color(dpg.mvThemeCol_WindowBg, (255, 255, 255), category=dpg.mvThemeCat_Core)
                                dpg.add_theme_color(dpg.mvThemeCol_TitleBg, (230, 230, 230), category=dpg.mvThemeCat_Core)
                                
                                dpg.add_theme_color(dpg.mvThemeCol_ChildBg, (255, 255, 255)) # Белый фон для дочерних окон (ВАЖНО для диалога!)
                                dpg.add_theme_color(dpg.mvThemeCol_PopupBg, (255, 255, 255)) # Фон всплывающих окон
                                dpg.add_theme_color(dpg.mvThemeCol_ModalWindowDimBg, (200, 200, 200, 150)) # Затемнение фона за модальным окном
                                                           
                                dpg.add_theme_color(dpg.mvThemeCol_Button, (0, 128, 255), category=dpg.mvThemeCat_Core)
                                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (180, 180, 180), category=dpg.mvThemeCat_Core)
                                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (255, 255, 255), category=dpg.mvThemeCat_Core)

                                dpg.add_theme_color(dpg.mvThemeCol_Text, (0, 0, 0), category=dpg.mvThemeCat_Core)

                                dpg.add_theme_color(dpg.mvThemeCol_TableHeaderBg, (255, 255, 255), category=dpg.mvThemeCat_Core)
                                dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (128, 128, 128), category=dpg.mvThemeCat_Core)

                # Применяем тему глобально
                dpg.bind_theme(light_theme)

        # Инициализация и установка темы
        #dpg.create_context()
        set_light_theme()

        #Определяем шрифт, который будет использоваться в программе
        with dpg.font_registry():
                with dpg.font(f"{path}/tahoma.ttf", 16, default_font=True, tag="Default font") as f:
                        dpg.add_font_range_hint(dpg.mvFontRangeHint_Cyrillic)
                        dpg.add_font_range(0x1F000, 0x1FFFF)
                        dpg.bind_font("Default font")                      

        #Определяем логику поведения файлового менеджера
        def file_dialog_callback(sender, app_data):
                global file_list
                global total_genomes, current_progress, progress_step

                if app_data is not None:
                        # Check if the user clicked "Ok"
                        if app_data['file_path_name'] is not None:
                                selected_paths = app_data['selections'] #list of paths
                                file_list = list(selected_paths.values())
                        else:
                                print("File selection cancelled.")

                # Количество геномов
                total_genomes = len(file_list)
                progress_step = 1.0 / total_genomes
                current_progress = 0.0	

        #Формирование основного интерфейса
        with dpg.window(label="Genome analysis", no_title_bar=True, tag="Primary Window", width=1000, height=400):
                with dpg.group(horizontal=True):
                        dpg.add_button(label="Выбрать файлы", callback=lambda: dpg.show_item("file_dialog"))
                        dpg.add_text("Формат выходного файла: ")
                        dpg.add_radio_button(
                                tag = "table_format_selector",
                                default_value = "tsv",                                
                                items = ["tsv", "json"],
                                horizontal = True
                        )
                        dpg.add_button(label="Анализ", callback=push)
                        dpg.add_loading_indicator(tag="loader", style=2, radius=1, color=(138, 43, 226))
                        dpg.hide_item("loader")                        
                        dpg.add_text("", tag="Analysis")

                dpg.add_progress_bar(tag="progress_bar", default_value=0.0, width=-1)
                dpg.add_text(tag="progress_text", default_value="Анализ: 0%", color=(0, 0, 0))		

                dpg.add_text("Валидация геномов", color=(0,0,255))
                with dpg.table(header_row=True, row_background=True, tag="quality_table", borders_outerH=True, borders_innerV=True, borders_innerH=True, borders_outerV=True):
                        dpg.add_table_column(label="Проба")
                        dpg.add_table_column(label="N ген")
                        dpg.add_table_column(label="G ген")
                        dpg.add_table_column(label="Валидация")

                dpg.add_text("Классификация геномов", color=(0,0,255))	
                with dpg.table(header_row=True, row_background=True, tag="class_table", borders_outerH=True, borders_innerV=True, borders_innerH=True, borders_outerV=True):
                        dpg.add_table_column(label="Проба")
                        dpg.add_table_column(label="Кластер по N гену")
                        dpg.add_table_column(label="Кластер по G гену")
                        dpg.add_table_column(label="Кластер по межд. класс.")

                with dpg.file_dialog(
                        width=450, height=300,
                        tag="file_dialog",
                        show=False,
                        callback=file_dialog_callback,
                        default_path=".",
                        file_count=0):
                        dpg.add_file_extension(".fasta", color=(0, 255, 255, 255), custom_text="[fasta]")
                        dpg.add_file_extension(".fa", color=(0, 255, 0, 255), custom_text="[fa]")              

        #Жизненный путь приложения		
        dpg.create_viewport(title=f'RabiesAnalyzerRus v.{version}', 
                            width=1000, height=400, small_icon='icon.png', large_icon='icon.png')
        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.start_dearpygui()
        dpg.destroy_context()

if __name__ == "__main__": #Проверяем, запускаемся ли мы в консоли или графике
        if len(sys.argv) <= 1:
                print("Запуск в GUI режиме")
                isGui = True
                run_gui()
        else:
                print("Запуск в CLI режиме")
                isGui = False	
                push()
