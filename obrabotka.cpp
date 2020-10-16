#include <iostream>
#include <iomanip> // для setw...
#include <string>  // подключаем строки
#include <fstream> // подключаем файлы
#include <stdlib.h>
#define _USE_MATH_DEFINES // для M_PI - числа пи
#include <math.h>
#include <complex> // для дифр поправки для коэф. затухания
#include <vector> // для функции автокорреляции (пока лишь)
#include <omp.h> // для быстрого подсчета автокорреляции

using namespace std; // используем стандартное пространство имен
//using namespace std::complex_literals;

unsigned __int64 findmax(double** a, unsigned __int64 size);         // поиск максимума

//void del1(double** a, unsigned int maxind, double delta);          // зануление вокруг максимума

void del(double** a, unsigned __int64 ind, unsigned __int64 deltai, unsigned __int64 maxind); // зануление вокруг максимума (по индексам)

double average(double* a, unsigned __int64 maxind);                  // среднее арифметическое одномерного массива

void sorta(double** a, unsigned __int64 maxind);                     // сортировка массива по возрастанию

void sortd(double** a, unsigned __int64 maxind);                     // сортировка массива по убыванию

void out(double** a, unsigned __int64 maxind);                       // вывод массива вида a[0][...], a[1][...]

double k_La(double f, double a, double c_L);                         // параметр k_La для расчета дифр. поправки; с_l - скорость

double z(double d, double a, unsigned __int64 n, unsigned __int64 m);                 // параметр z для расчета дифр. поправки; n,m  - номера пиков

double c_difr(double c_L, double k_La, double z);                     // дифракционная поправка к скорости

complex <double> difrhelp(double k_La, double d, double a, unsigned __int64 n);       // вспомогательная ф-я для нахождения дифр. поправки к коэф. затухания

double difr(double k_La, double d, double a, unsigned __int64 n, unsigned __int64 m); // дифр. поправка для коэф. затухания

double iterative_avg(double** a, unsigned __int64 length, unsigned int k, double eps, unsigned int n); //функция для подсчета "среднего" уровня шума
    
double autocorrelation(double** a, unsigned __int64 j, unsigned __int64 n);           // функция автокорреляции сигнала

unsigned __int64 countlines(string filename);                                         // подсчет строк в файле

void fill_u(double** array, unsigned __int64 maxind, string filename, string regime); // чтение файла и заполнение массива

void initialize_u(double**& array, unsigned __int64 maxindex);         //инициализация массивов array[0][...], array[1][...]

void autocorrelation_method(double** array, unsigned __int64 maxindex);               // вычисляем функцию автокорреляции сигнала, записываем в файл

unsigned __int64 delta_in_indices(double delta, double duration, unsigned __int64 maxindex); // преобразование "ширины" импульса в ширину "по индексам", т.е. время delta занимает столько-то индексов

unsigned __int64 deltai_modification(unsigned __int64 deltai, double frequency);      // зависимость deltai (ширины импульса) от неких параметров, пока что только от частоты

double k_modification(double k, double f); // зависимость k (что-то типа отношения сигнал/шум) от внешних параметров

void filter_noise(double** array, unsigned __int64 maxindex, double k, double avg);   // "фильтрация" массива амплитуд от "шума"

void fill_u_max(double** array_fill_from, unsigned __int64 maxindex_from, double** array_to_be_filled, unsigned __int64 maxindex_to_be_filled, unsigned __int64 deltai);
// заполнение максимумами промежуточного массива u_max и подсчет ненулевых элементов в нем, чтобы сформировать окончательный массив

void check_peak_distance(double** array, unsigned __int64 maxindex, double peak_distance); // проверка расстояний между пиками

void calculate_speed(double* array_to_be_calculated, double** array_to_calculate_from, unsigned __int64 maxindex, double frequency, double specimen_length, double electrode_diameter);
// вычисление скорости с учетом дифракционной поправки

void calculate_attenuation(double* speed_array, double* array_to_be_calculated, double** array_to_calculate_from, unsigned __int64 maxindex, double frequency, double specimen_length, double electrode_diameter);
// расчет коэф. затухания с учетом дифр. поправки

void file_out(double** final_array, double* speed_array, double* attenuation_array, unsigned __int64 maxindex, double specimen_length, double frequency, double electrode_diameter, string filename, string regime);
// вывод в файл

void compute(string filename, string filename2, string regime, double d, double delta, double f, double a); // почти главная функция

double norma(double** array, unsigned __int64 maxindex); // норма функции сигнала u(t) для нормирования ф-и автокорреляции

unsigned __int64 find_corresponding_amplitude(double** array_where_to_find, unsigned __int64 maxindex, double corresponding_time); // поиск амплитуды для соотв. интервала времени в автокорреляционном методе; возвращает индекс

//unsigned __int64 Q = 0;

double firstmax[2] = { 0,0 }; // костыль; первый максимум в исходных данных; нужен для расчета коэф. заутхания (конкретно для поиска соотв. амплитуды)
// firstmax[0] - время, firstmax[1] - амплитуда


int main() {
    
    
    //string regime = "autocorrelation";
    //------------------------------------------------------------------------------------------------------------

    //string filename = "autocorrelation.txt";
    //string filename = ".\\waveforms\\C1Quartz_19,948 mm_5 MHz_00000.txt";
    string filename = ".\\waveforms\\F1_D16_15,02 mm_5 MHz_00000.txt";    
    //string filename = ".\\waveforms\\C1_orgsteklo-20,454-2,5MHz_00000.txt";    
    //string filename = ".\\waveforms\\C1_orgsteklo kvadrat-19,589-5MHz_00000.txt";

    cout << "file name: " << filename << endl;

    // внешние данные
    //------------------------------------------------------------------------------------------------------------

    double d; // толщина (длина) образца
    cout << "enter the d, meters:" << endl << "d=";
    //cin >> d;
    //d = 0.019948;
    d = 0.01502;
    //d = 0.020454;
    cout << endl;

    double delta; // длина радиоимпульса по времени
    delta = 2*310e-9; // 5 периодов
    cout << endl;

    double f; // частота, Гц
    cout << "enter the frequency, Hz" << endl << "f=";
    cin >> f;
    //f = 2.5e6;
    cout << endl;

    double a; // радиус электрода в метрах
    a = 0.01;
    cout << endl;

    //------------------------------------------------------------------------------------------------------------
    
    compute(filename, filename, "findmax", d, delta, f, a); // ищем максимумы и вдобавок вычисляем ф-ю автокорреляции
    compute("autocorrelation.txt", filename, "autocorrelation", d, delta, f, a); // используя вычисленную ф-ю автокорреляции, ищем максимумы её и получаем скорость

   //cout << "Q =" << Q;
    

    system("pause");

    return 0;
}


unsigned __int64 findmax(double** a, unsigned __int64 size) // ищем глобальный максимум амплитуды
{
    double max = a[1][0];
    unsigned __int64 maxindex = 0;
    unsigned __int64 i; //номер макс элемента
    for (i = 1; i < size; i++) {
        if (a[1][i] == 0) continue; // есть ли смысл?
        if (a[1][i] > max) {
            max = a[1][i];
            maxindex = i;
        }
    }
    return maxindex;
}

void del1(double** a, unsigned __int64 maxind, double delta) // зануляем все элементы вокруг максимума в радиусе delta (по времени)
{
    unsigned __int64 i=0;
    while (a[0][i] > (a[0][maxind] - delta) && a[0][i] > (a[0][maxind] + delta)) {
        a[1][i] = 0;
        i++;
    }
}

void del(double** a, unsigned __int64 ind, unsigned __int64 deltai, unsigned __int64 maxind) // другой вариант (лучше, вроде) зануления всех элементов вокруг выбранного элемента в радиусе deltai (по индексам)
{
    unsigned __int64 upper_bound = ind + deltai;
    unsigned __int64 lower_bound = ind - deltai; // при unsigned int возникает ошибка, если получается отрицательное (ind-deltai) при использовании эквивалентного условия

    if (ind < deltai) // <=> (ind - deltai < 0)
        lower_bound = 0;

    if (upper_bound > maxind)
        upper_bound = maxind;

    //for (unsigned __int64 i = ind - deltai; i < ind + deltai; i++)
    for (unsigned __int64 i = lower_bound; i < upper_bound; i++)
        a[1][i] = 0;
}

double average(double* a, unsigned __int64 maxind) // подсчет среднего значения одномерного массива
{
    double avgvalue = 0;
    for (unsigned __int64 i = 0; i < maxind; i++)
        avgvalue += a[i];

    return avgvalue / maxind;
}

void sorta(double** a, unsigned __int64 maxind) { // сортировка массива пузырьком по возрастанию
    double temp0, temp1; //временные переменные для обмена, соотв. a[0][...] и a[1][...]
    for (unsigned __int64 i = 0; i < maxind - 1; i++) {
        for (unsigned __int64 j = 0; j < maxind - i - 1; j++) {
            if (a[0][j] > a[0][j + 1]) {
                temp0 = a[0][j];
                temp1 = a[1][j];

                a[0][j] = a[0][j + 1];
                a[1][j] = a[1][j + 1];

                a[0][j + 1] = temp0;
                a[1][j + 1] = temp1;
            }
        }
    }
}

void sortd(double** a, unsigned __int64 maxind) { // сортировка массива пузырьком по убыванию
    double temp0, temp1; //временные переменные для обмена, соотв. a[0][...] и a[1][...]
    for (unsigned __int64 i = 0; i < maxind - 1; i++) {
        for (unsigned __int64 j = 0; j < maxind - i - 1; j++) {
            if (a[0][j] < a[0][j + 1]) {
                temp0 = a[0][j];
                temp1 = a[1][j];

                a[0][j] = a[0][j + 1];
                a[1][j] = a[1][j + 1];

                a[0][j + 1] = temp0;
                a[1][j + 1] = temp1;
            }
        }
    }
}

void out(double** a, unsigned __int64 maxind) { // вывод массива вида a[0][...], a[1][...]    
    for (unsigned __int64 i = 0; i < maxind; i++) {
        cout << setw(6) << scientific << setprecision(6) << a[0][i] << ", " << resetiosflags(ios_base::floatfield);
        cout << setw(6) << setprecision(6) << a[1][i] << endl;
    }
    cout << endl << endl;
    cout << resetiosflags(ios_base::floatfield);
}

double k_La(double f, double a, double c_L) {  // параметр k_La для расчета дифр. поправки; с_l - скорость
    return 2 * M_PI * f * a / c_L;
} 

double z(double d, double a, unsigned __int64 n, unsigned __int64 m) { // параметр z для расчета дифр. поправки; n,m  - номера пиков
    return d / a * sqrt((2 * n - 1) * (2 * m - 1));
}

double c_difr(double c_L, double k_La, double z) { // дифракционная поправка к скорости
    return -c_L * (5.2 / (pow(z, 3 / 2) * k_La * k_La) - 7e-4 * (k_La * k_La - 2200) / (k_La * k_La + 13 * z * z));
}

complex<double> difrhelp(double k_La, double d, double a, unsigned __int64 n) { // вспомогательная ф-я для нахождения дифр. поправки к коэф. затухания

    //double xi = k_La / (2 * a) * floor(sqrt((2 * n - 1) * (2 * n - 1) * d * d + 4 * a * a) - (2 * n - 1) * d); // параметр кси для расчета
    //cout << "difrhelp: k=" << k_La / a << endl;
    double xi = k_La / (2 * a) * (sqrt((2 * n - 1) * (2 * n - 1) * d * d + 4 * a * a) - (2 * n - 1) * d); // без floor
    cout << "difrhelp: xi=" << xi << endl;
    return 1.0 - ((1 - xi * xi / (2 * k_La * k_La)) * _j0(xi) + complex<double>(0, 1) * (1 - xi * xi / (2 * k_La * k_La) + xi / (k_La * k_La)) * _j1(xi)) * exp(-complex<double>(0, 1) * xi);
    }

double difr(double k_La, double d, double a, unsigned __int64 n, unsigned __int64 m) { // дифр. поправка для коэф. затухания
    return 20 * log10(abs(difrhelp(k_La, d, a, n) / difrhelp(k_La, d, a, m)));
    }

double iterative_avg(double** a, unsigned __int64 length, unsigned int k, double eps, unsigned int n) { //функция для подсчета "среднего" уровня шума
    // a - массив (одномерный) для обработки
    // length - длина массива a
    // k - коэф. для сравнения среднего и текущего значений; используется для отфильтровывания больших значений, т.е. не шума
    // eps - число для оценки изменений нового и старого средних; "точность"
    // n - максимальное число итераций

    double avg = 0;
    for (unsigned __int64 i = 0; i < length; i++)
        avg += abs(a[1][i]);
    avg = avg / length;

    cout << "iterative avg: avg=" << avg << endl;

    double newavg = 0;
    
    for (unsigned __int64 i = 0; i < length; i++) { // считаем newavg для использования в условии ниже
        if (abs(a[1][i]) > avg * k)
            newavg += avg * k;
        else newavg += abs(a[1][i]);
    }
    newavg = newavg / length;

    cout << "iterative avg: newavg=" << newavg << endl;
    
    double ratio = avg / newavg; //для проверки условия

    cout << "iterative avg: initial avg/newavg=" << avg / newavg << endl;

    for (unsigned int i = 0; ((i < n) && (ratio > 1 + eps)); i++) {

        cout << "iterative avg: iteration #" << i << endl;

        newavg = 0;

        for (unsigned __int64 i = 0; i < length; i++) { // на основе ранее вычисленного среднего вычисляем новое "модифицированное" среднее
            if (abs(a[1][i]) > avg * k)
                newavg += avg * k;
            else newavg += abs(a[1][i]);
        }
        newavg = newavg / length;

        cout << "iterative avg: newavg=" << newavg << endl;

        ratio = avg / newavg;

        cout << "iterative avg: avg/newavg=" << ratio << endl;

        avg = newavg;   

        
    }
    return newavg;
}

double autocorrelation(double** a, unsigned __int64 j, unsigned __int64 n) { // функция автокорреляции сигнала
    // a - массив {t,A}
    // j - сдвиг по времени (индексу)
    // n - число элементов в дискретном сигнале

    //double R = 0;
    /*
    for (unsigned __int64 i = 0; i < n-j; i++) { // чушь какая-то написана здесь в R+=
        //cout << "a=" << a[1][i];
        R += a[1][i] * a[1][i + j] * (a[0][i + j] + a[0][i]); 
        //R += a[1][i] * a[1][i + j];
    }    

   R = R / (a[0][n] - a[0][0]); 
   cout << "R=" << R << endl;*/

    double Rt = 0;

    /*for (unsigned __int64 i = 0; i < n - j - 1; i++) { // метод левых прямоугольников
        R += a[1][i] * a[1][i + j] * (a[0][i + 1] - a[0][i]);        
    }
    R = R / (a[0][n] - a[0][0]);

    cout << "R=" << R << endl;*/
    
    for (unsigned __int64 i = 0; i < n - j - 1; i++) { // метод трапеций
        Rt += (a[1][i] * a[1][i + j] + a[1][i + 1] * a[1][i + j + 1]) / 2 * (a[0][i + 1] - a[0][i]);
        //Q++;
    }
    Rt = Rt / (a[0][n] - a[0][0]);
    
    // значение по методу трапеций, опять же, очень близко к значению по методу прямоугольников

    //cout << "Rt=" << Rt << endl;

    //cout << "R/Rt=" << R/Rt << endl;

    return Rt;
}

unsigned __int64 countlines(string filename) { //подсчет строк в файле
    //открываем файл и считаем количество строчек
    ifstream file(filename);
    if (!file) { // проверка открытия файла
        cout << "error opening file" << endl;
        return -1; 
    }
    unsigned __int64 counter = 0;
    while (file.ignore((numeric_limits<streamsize>::max)(), '\n')) // пропускаем и считаем строки
        counter++;
    file.close();
    //cout << "counter= " << counter << endl;
    //посчитали строчки
    return counter;
}

void fill_u(double** array, unsigned __int64 maxind, string filename, string regime) { // чтение файла и заполнение массива

    //считываем значения и записываем в массив
    ifstream file1(filename); // файл, из которого читаем
    string s; // сюда будем класть считанные строки
    if (!file1)
        cout << "error reading file" << endl;


    if (regime == "findmax") { //игнорируем  первые 5 строк, если режим поиска максимумов
        file1.ignore((numeric_limits<streamsize>::max)(), '\n');
        file1.ignore((numeric_limits<streamsize>::max)(), '\n');
        file1.ignore((numeric_limits<streamsize>::max)(), '\n');
        file1.ignore((numeric_limits<streamsize>::max)(), '\n');
        file1.ignore((numeric_limits<streamsize>::max)(), '\n');
    }

    //считываем cтроки с файла и записываем в массив

//#pragma omp parallel for  schedule(static)
    for (unsigned __int64 j = 0; j < maxind; j++) { // -1, т.к. последняя строка - пустая /эээ

        getline(file1, s);
        auto pos = s.find(",");
        if (pos != string::npos)
        {
            array[0][j] = stod(s.substr(0, pos)); //string to double
            array[1][j] = stod(s.substr(pos + 1));
        }
    }
    file1.close(); // закрываем файл
    //считали, заполнили
}

void initialize_u(double**& array, unsigned __int64 maxindex) { //инициализация массивов "u"      
    array = new double* [2];

    array[0] = new double[maxindex];
    array[1] = new double[maxindex];
}

void autocorrelation_method(double** array, unsigned __int64 maxindex) { // вычисляем функцию автокорреляции сигнала, записываем в файл
    
    double** acor = new double* [2];
    acor[0] = new double[maxindex];
    acor[1] = new double[maxindex];

    double u_norma = norma(array, maxindex);

    cout << "u_norma=" << u_norma << endl;

    //cout << "autocor=" << autocorrelation(u, 0, maxindex);

#pragma omp parallel for schedule (dynamic, 100) //num_threads(8)    //dynamic, 1000
    //for ( unsigned __int64 i = 0; i < counter - 1; i++) {
    // записываем в массив acor результаты ф-и автокорреляции по аналогии с массивом u: acor[0][...] - время,  acor[1][...] - значение ф-и автокор.
    for (__int64 i = 0; i < static_cast<__int64>(maxindex); i++) {     
        acor[0][i] = array[0][i] - array[0][0];        
        acor[1][i] = autocorrelation(array, i, maxindex) / u_norma; // вычисляем функцию автокорреляции R(tau)=1/T int_0^T  u(t)*u(t+tau) dt
        
        //cout << "acor=" << acor[1][i] << endl;
    }

    // выводим в файл

    ofstream corfile("autocorrelation.txt"); // открываем файл для записи

    for (unsigned __int64 i = 0; i < maxindex; i++)
        corfile << setprecision(12) << acor[0][i] << "," << acor[1][i] << endl; // записываем
    corfile.close();
}

double norma(double** array, unsigned __int64 maxindex) { // норма функции сигнала u(t) для нормирования ф-и автокорреляции
    // вычисляем интеграл N=  1/{t=t(maxindex)} * int_{t=0}^{t=t(maxindex)} u(t)^2 dt методом трапеций

    double N = 0;
    for (unsigned __int64 i = 0; i < maxindex - 1; i++) {
        N += (array[1][i] * array[1][i] + array[1][i + 1] * array[1][i + 1]) / 2 * (array[0][i + 1] - array[0][i]);
    }
   
    return N / (array[0][maxindex] - array[0][0]);
}

unsigned __int64 delta_in_indices(double delta, double duration, unsigned __int64 maxindex) { // преобразование "ширины" импульса в ширину "по индексам", т.е. время delta занимает столько-то индексов
    // duration - длительность осциллограммы в секундах
    // maxindex - длина массива

    return static_cast<unsigned __int64> ((delta * maxindex) / duration);
}

unsigned __int64 deltai_modification(unsigned __int64 deltai, double frequency) { // зависимость deltai (ширины импульса) от неких параметров, пока что только от частоты
    return static_cast<unsigned __int64> (ceil(deltai * 15 * pow(frequency * 1e-6, -0.8))); // чем меньше частота, тем "шире" пик, и для компенсации этого увеличиваем k на низких частотах
}

double k_modification(double k, double f) { // зависимость k (что-то типа отношения сигнал/шум) от внешних параметров, пока что от частоты
    return k *= -0.02 * f * f * 1e-12 + 1.15 * f * 1e-6 - 1; // ~ чем больше частота, тем "выше" пик    
}

void filter_noise(double** array, unsigned __int64 maxindex, double k, double avg) { // "фильтрация" массива амплитуд от "шума" (небольших амплитуд) на основе коэф. k и "среднего" значения        
    for (unsigned __int64 i = 0; i < maxindex; i++) {
        //u[1][i] = abs(u[1][i]); // можно применить из-за несимметричности (иногда) осциллограммы относительно оси абсцисс
        array[1][i] -= avg * k;
        if (array[1][i] < 0) // если за вычетом avg*k получилось отриц. значение, то зануляем его, и в дальнейшем интересовать нулевые значения не будут
            array[1][i] = 0; 
    }
}

void fill_u_max(double** array_fill_from, unsigned __int64 maxindex_from, double** array_to_be_filled, unsigned __int64 maxindex_to_be_filled, unsigned __int64 deltai) {
    // заполнение максимумами промежуточного массива u_max и (убрано) подсчет ненулевых элементов в нем, чтобы сформировать окончательный массив

    // заполняем массив с максимумами:
    // 1) находим глобальный max
    // 2) записываем его в соотв. массив - array_to_be_filled
    // 3) зануляем элементы вокруг максимума
    // (УБРАНО)   // 4) считаем количество (m) ненулевых элементов из array_to_be_filled, чтобы потом сфорировать окончательный массив
    //unsigned __int64 m = 0; // длина будущего массива u_m "окончательного"
    for (unsigned __int64 i = 0; i < maxindex_to_be_filled; i++) {
        unsigned __int64 max = findmax(array_fill_from, maxindex_from);
        //cout << "max=" << max << endl;
        array_to_be_filled[0][i] = array_fill_from[0][max];
        array_to_be_filled[1][i] = array_fill_from[1][max]; // значение за вычетом k*avg    
        //if (array_to_be_filled[1][i] > 0) m++;
        del(array_fill_from, max, deltai, maxindex_from);
    }    
}

void check_peak_distance(double** array, unsigned __int64 maxindex, double peak_distance) { // проверка расстояний между пиками
    /*for (unsigned __int64 i = maxindex; i > 0; i--) { // обходим массив и проверяем расстояние между пиками:
                                                      // если расстояние меньше 0,9 * peak_distance , то пик отсеивается;
                                                      //0,9 просто так
        if (array[1][i] == 0) continue;
        if (array[0][i] - array[0][i - 1] < peak_distance * 0.9) {
            array[1][i - 1] = 0;
        }
    }*/

    unsigned __int64 global_max = findmax(array, maxindex); // ищем глобальный максимум, чтобы от него отсчитывать расстояние потом

    for (unsigned __int64 i = global_max; i < maxindex; i++) { // вперед от максимума
        if (array[1][i] == 0) continue;
        if (array[0][i + 1] - array[0][i] < peak_distance * 0.9)
            array[1][i + 1] = 0;
    }

    for (unsigned __int64 i = global_max; i > 0; i--) { // назад от максимума
        if (array[1][i] == 0) continue;
        if (array[0][i] - array[0][i - 1] < peak_distance * 0.9)
            array[1][i - 1] = 0;
    }
}

/*void check_peak_distance2(double** array, unsigned __int64 maxindex) { // проверка расстояний между пиками

    double* distance_array; // массив попарных расстояний
    distance_array = new double[maxindex * (maxindex - 1) / 2]; // размер - число пар без повторений из maxindex элементов

    unsigned __int64* hits_array; // массив "попаданий" для каждой точки
    hits_array = new unsigned __int64[maxindex];

    for (unsigned __int64 i=0; i<maxindex; i++) //строим все (без повторений) пары расстояний между точками
        for (unsigned __int64 j = 0; j < maxindex; j++) 
            if (i > j) { // "выше главной диагонали"
                distance_array[i + j] = array[1][i] - array[1][j];
            }
    
    for (unsigned __int64 i=0; i<maxindex; i++) // каждая точка из массива
        for (unsigned __int64 j = 0; j < (maxindex * (maxindex - 1) / 2); j++) { // каждое расстояние
            array[0][i]+distance_array
        }

}*/

void calculate_speed(double* array_to_be_calculated, double** array_to_calculate_from, unsigned __int64 maxindex, double frequency, double specimen_length, double electrode_diameter) { // расчет скорости с учетом дифракционной поправки
    // array_to_be_calculated - массив, значения которого вычисляются (массив скоростей)
    // array_to_calculate_from - массив, на основе значений которого вычисляется массив выше
    // maxindex - длина массива выше
    // frequency - частота, нужна для вычисления дифракционной поправки
    // specimen_length - длина (толщина) образца
    // electrode_diameter - диаметр электрода, нужен для выч. дифр. поправки

    for (unsigned __int64 i = 0; i < maxindex - 1; i++) {
        array_to_be_calculated[i] = 2 * specimen_length / (array_to_calculate_from[0][i + 1] - array_to_calculate_from[0][i]);
        //cout << endl << "deltaC_dif=" << c_dif(array_to_be_calculated[i], k_La(frequency, electrode_diameter, array_to_be_calculated[i]), z(specimen_length, electrode_diameter, i + 1, i + 2)); // для теста
        array_to_be_calculated[i] += c_difr(array_to_be_calculated[i], k_La(frequency, electrode_diameter, array_to_be_calculated[i]), z(specimen_length, electrode_diameter, i + 1, i + 2)); // с дифракционной поправкой
        cout << "v[" << i << "] = " << array_to_be_calculated[i] << endl;
    }
}

void calculate_attenuation(double* speed_array, double* array_to_be_calculated, double** array_to_calculate_from, unsigned __int64 maxindex, double frequency, double specimen_length, double electrode_diameter) { // расчет коэф. затухания с учетом дифр. поправки
    
    // speed_array - массив скорости (для k_La)
    // array_to_be_calculated - массив, значения которого вычисляются (массив коэф. затухания)
    // array_to_calculate_from - массив, на основе значений которого вычисляется массив выше
    // maxindex - длина массива выше
    // frequency - частота, нужна для вычисления дифракционной поправки
    // specimen_length - длина (толщина) образца (d)
    // electrode_diameter - диаметр электрода, нужен для выч. дифр. поправки

    for (unsigned __int64 i = 0; i < maxindex - 1; i++) {
        double Adifr = difr(k_La(frequency, electrode_diameter, speed_array[i]), specimen_length, electrode_diameter, i + 1, i + 2); // дифр поправка
        cout << "calculate_attenuation: difr. correction (dB) = " << Adifr << endl;
        array_to_be_calculated[i] = (20 * log10(array_to_calculate_from[1][i] / array_to_calculate_from[1][i + 1]) /*- Adifr*/) / (2 * specimen_length);
        cout << "alpha[" << i << "] = " << array_to_be_calculated[i] << endl;
    }
}

void file_out(double** final_array, double* speed_array, double* attenuation_array, unsigned __int64 maxindex, double specimen_length, double frequency, double electrode_diameter, string filename, string regime) { // вывод в файл

    // final_array - массив с найденными максимумами
    // speed_arra - м. скоростей
    // attenuation_array - м. с коэф. затухания
    // maxindex - длина final_array (у других на 1 меньше)

    string outfile_name;

    if (regime=="findmax")
        outfile_name = "results_findmax.txt";

    if (regime == "autocorrelation")
        outfile_name = "results_autocorrelation.txt";

    ofstream outfile(outfile_name);

    outfile << filename << endl << "d=" << specimen_length << " m" << endl << "f=" << frequency << " Hz" <<endl << "a=" << electrode_diameter << " m" << endl;

    outfile << endl;

    for (unsigned __int64 i = 0; i < maxindex; i++) {
        outfile << "max #" << i << " = " << scientific << setprecision(16) << final_array[0][i] << ", " << resetiosflags(ios_base::floatfield);
        outfile << setprecision(14) << final_array[1][i] << endl << resetiosflags(ios_base::floatfield);
    }

    /*if (regime == "findmax") {
       outfile << endl;

        for (unsigned __int64 i = 0; i < maxindex - 1; i++)
            outfile << setprecision(6) << "alpha[" << i << "] = " << attenuation_array[i] << " dB/m" << endl;
    }*/

    for (unsigned __int64 i = 0; i < maxindex - 1; i++)
        outfile << setprecision(6) << "alpha[" << i << "] = " << attenuation_array[i] << " dB/m" << endl;

    outfile << endl << endl;

    for (unsigned __int64 i = 0; i < maxindex - 1; i++)
        outfile << setprecision(6) << "v[" << i << "] = " << speed_array[i] << " m/s" << endl;

    outfile.close();
}

unsigned __int64 find_corresponding_amplitude(double** array_where_to_find, unsigned __int64 maxindex, double corresponding_time) { // поиск амплитуды для соотв. интервала времени в автокорреляционном методе; возвращает индекс
    // array_where_to_find - массив, в котором искать соотв амплитуду
    // maxindex - длина массива выше
    // corresponding_time - временной интервал, для которого надо найти амплитуду

    // берем массив для поиска, и, сравнивая время, ищем индекс массива, соотв. времени
    
    unsigned __int64 i=0; // счётчик и возвращаемый индекс

    while (array_where_to_find[0][i] < corresponding_time)
        i++;
    cout << "find cor ampl: i=" << i << ", time= " << array_where_to_find[0][i]<< " ampl=" << array_where_to_find[1][i] << endl;
    
    return i;
}

void compute(string filename, string filename2, string regime, double d, double delta, double f, double a) { // вычисляем всю фигню
    // filename2 - костыль - дополнительный файл (и массив), содержащий исходные данные; нужен для расчета коэф. затухания в автокорреляционном методе

    //------------------------------------------------------------------------------------------------------------

    
    cout << "*****  " << regime << "  *****" << endl << endl;

    double** u; //массив для хранения времени сигнала и его величины: u[0] - массив времён, u[1] - массив амплитуд


    unsigned __int64 counter = countlines(filename); // открываем файл и считаем количество строк
    cout << "countlines: counter = " << counter << endl;


    // counter будет в дальнейшем использован для чтения из файла в цикле; первые 5 строчек - ненужная (пока) информация   



    if (regime == "findmax") // если поиск максимумов, то уменьшаем длину массива на 5, т.к. 5 строк в файле осциллограммы содержат техническую информацию
        counter -= 5;

    //counter--; // еще в конце файла последняя строка - пустая

    initialize_u(u, counter - 1); // выделение памяти под массивы u[0][...], u[1][...]

    fill_u(u, counter - 1, filename, regime);

    //------------------------------------------------------------------------------------------------------------
    // лепим костыль по аналогии с созданием массива 'u'
    double** u2;
    
    unsigned __int64 counter2 = countlines(filename2);
    counter2 -= 5;
    cout << "counter2=" << counter2 << endl;
    initialize_u(u2, counter2 - 1);
    //if (regime=="autocorreltation")
        fill_u(u2, counter2 - 1, filename2, "findmax");

    //out(u2, counter2);
    
    //слепили
    //------------------------------------------------------------------------------------------------------------


    if (regime == "findmax") // если режим поиска максимумов, то вычисляем ф-ю автокорреляции для использования потом, и ищем максимумы
        autocorrelation_method(u, counter - 1);

    cout << "counter=" << counter << endl; //должно быть больше на 1 чем число, записанное в заголовке файла с осциллограммой

    //выводим первый и последний элементы
    cout << "u[0]=" << u[0][0] << ", " << u[1][0] << endl;
    cout << "index of the last element is " << counter - 2 << ", u[last]=" << u[0][counter - 2] << ", " << u[1][counter - 2] << endl;


    //out(u, counter - 1); //вывод всех элементов массива

    //------------------------------------------------------------------------------------------------------------



    cout << "counter-2=" << counter - 2 << "    u[0][counter-2] = " << u[0][counter - 2] << endl;

    double tlen = u[0][counter - 2] - u[0][0]; //длительность осциллограммы по времени

    cout << endl << "tlen=" << tlen << endl;

    cout << "delta=" << delta << " seconds" << endl;




    unsigned __int64 deltai = delta_in_indices(delta, tlen, counter - 1); // переводим длительность импульса по времени delta в длину "по индексам"

    cout << "deltai=" << deltai << " indices" << endl;

    deltai = deltai_modification(deltai, f); // модифицируем deltai  в зависиомости от частоты: меньше частота - шире пик
    cout << "modified deltai=" << deltai << " indices" << endl;

    delta = deltai * tlen / (counter - 1);
    cout << "(~ in seconds) modified delta=" << delta << endl;


    //------------------------------------------------------------------------------------------------------------


    // вычисляем среднюю амплитуду (по абсолютному значению)
    // считаем, не используя ф-ю average, т.к. интересует среднее по абс. значениям
    double avg = 0;
    for (unsigned __int64 i = 0; i < counter - 1; i++)
        avg += fabs(u[1][i]); //вроде правильно, fabs принимает и выдает double
    avg = avg / (counter - 1);

    cout << "avg=" << avg << endl;
    // вычислили среднее, вывели на экран
    // среднее значение, вычисленное по методу трапеций, не отличается значительно от того, что выше


    // вычислияем "среднее значение" так:
    //      1) вычисляем среднее как выше
    //      2) используя полученное значение, отсеиваем большие по сравнению с этим средним значения
    //      3) получаем новое среднее значение
    //      4) шаг 2
    //      5) и так далее, пока не вступит ограничение на число итераций или на относительную разницу между "новым" и "старым" средними значениями
    // т.о. образуется сходящаяся (?) последовательность значений
    double itavg = iterative_avg(u, counter - 1, 5, 0.05, 10);
    // аргументы: u - массив для обработки
    // counter - 1 - длина массива
    // 5 - коэф. для отсеивания больших значений: отсеиваются значения, в 5 раз больше среднего
    // 0.05 - ограничение на относ. разницу нового и старого средних знач.: если отн. разница меньше 0,05, то итерации прекращаются
    // 10 - ограничение на число итерация: выход из цикла после 10 итерации


    cout << "iterative avg=" << itavg << endl << "avg/itavg=" << avg / itavg << endl;

    avg = itavg;

    //------------------------------------------------------------------------------------------------------------

    double k = 1; //как бы отношение сигнал/шум
    cout << endl << "initial k=" << k << endl;

    k = k_modification(k, f); // модифицируем k в зависимости от частоты
    k *= 0.25;
    cout << "k after modification: k=" << k << endl;

    filter_noise(u, counter - 1, k, avg); // отнимаем от всех амлитуд k средних значений, как бы фильтруя сигнал от шума


    //вывод ненулевых элементов массива
    /*
    for (int i = 0; i < counter - 1; i++) {
        if (u[1][i] == 0) continue;
        cout << "u[" << i << "]=" << u[0][i] << ", " << u[1][i] << endl;

    }
    */

    const unsigned __int64 l = 20; // длина массива с максимумами
    double** u_max; // массив с максимумами (не окончательный)
    initialize_u(u_max, l); // выделение памяти под массив u_max



    //------------------------------------------------------------------------------------------------------------

    unsigned __int64 gmax = findmax(u, counter - 1); //глобальный максимум
    cout << "gmax=" << gmax << ", " << "glob max=" << u[0][gmax] << ", " << u[1][gmax] + k * avg << endl;



    fill_u_max(u, counter - 1, u_max, l, deltai); // l - длина массива с максимумами
    // заполняем массив с максимумами:
    // 1) находим глобальный max
    // 2) записываем его в соотв. массив
    // 3) зануляем элементы вокруг максимума


    unsigned __int64 m = 0; // число ненулевых элементов в u_max
    for (unsigned __int64 i = 0; i < l; i++) // подсчитываем число ненулевых элементов в u_max
        if (u_max[1][i] > 0)
            m++;

    cout << endl << "m=" << m << endl;

    cout << "u_max=" << endl;
    out(u_max, l); // значение амплитуды выводится за вычетом k*avg
    cout << endl;


    //------------------------------------------------------------------------------------------------------------    
    // сортируем пузырьком u_max по времени по (!) убыванию, чтобы узнать расстояние между пиками и затем избавиться от первого пика, если он не годится

    sortd(u_max, m);

    cout << "descending u_max=" << endl;
    out(u_max, l);




    //вычисляем среднее расстояние (по времени) между максимумами, чтобы, если случайно 2 макс. "плохие", то это не нарушило бы ничего
    /*r = 0;
    for (int i = 0; i < m - 1; i++) {
        r += u_max[0][i] - u_max[0][i + 1];
    }
    r = r / (m-1);
    cout << endl << "r=" << r << endl;  */

    double* peak_distance = new double[m - 1];
    for (int i = 0; i < m - 1; i++) {
        peak_distance[i] = u_max[0][i] - u_max[0][i + 1];
    }
    double avg_peak_distance = average(peak_distance, m - 1);
    cout << "average peak distance = " << avg_peak_distance << endl;


    //------------------------------------------------------------------------------------------------------------

    cout << endl;
    sorta(u_max, m); // сортируем пузырьком u_max по времени по возрастанию, чтобы  потом избавиться от ненужного пика

    cout << "u_max before del=" << endl;
    out(u_max, l);

    //превратить в ф-ю
    /*for (unsigned __int64 i = m; i > 0; i--) { // обходим массив u_max и проверяем расстояние между пиками:
                                               // если расстояние меньше 0,9*r , то пик отсеивается; 0,9 просто так
        if (u_max[1][i] == 0) continue;
        if (u_max[0][i] - u_max[0][i - 1] < r * 0.9) {
            u_max[1][i - 1] = 0;
        }
    }*/

    check_peak_distance(u_max, m, avg_peak_distance);


    cout << "u_max after del=" << endl;
    out(u_max, l);

    m = 0; //пересчитываем m - количество ненулевых элементов из u_max, чтобы потом сфорировать окончательный массив
    for (unsigned __int64 i = 0; i < l; i++)
        if (u_max[1][i] > 0) m++;
    cout << endl << "m=" << m << endl;


    //------------------------------------------------------------------------------------------------------------

    cout << "ascending u_max=" << endl;
    out(u_max, l);

    double** u_m; // окончательный массив
    initialize_u(u_m, m);

    //записываем элементы из u_max в u_m  
    unsigned __int64 j = 0;
    for (unsigned __int64 i = 0; i < l; i++) {
        if (u_max[1][i] != 0) {
            u_m[0][j] = u_max[0][i];
            if (regime =="findmax")
                u_m[1][j] = u_max[1][i] + avg * k; // прибавляем k*avg, потому что раньше отнимали
            if (regime == "autocorrelation")
                u_m[1][j] = u2[1][find_corresponding_amplitude(u2, counter2, firstmax[0] + u_max[0][i])]; //+ avg * k; // firstmax - костыль; по идее - первый максимум
            j++;
        }
    }
    // записали

    // выводим u_m
    cout << endl << "u_m=" << endl;
    out(u_m, m);
    cout << endl << endl;

    //cout << "u_m00=" << u_m[0][0] << ", " << u_m[1][0] << endl;


    //------------------------------------------------------------------------------------------------------------
    /*double* hi; // логарифмический декремент затухания
    hi = new double[m];
    double avghi = 0; //среднее значение декремента

    for (unsigned __int64 i = 0; i < m - 1; i++) {
        hi[i] = log( u_m[1][i] / u_m[1][i + 1] );
        cout << "hi[" << i << "] = " << hi[i] << endl;
    }


    avghi = average(hi, m - 1);
    cout << endl << "avghi=" << avghi << endl << endl;*/
    // посчитали и вывели декремент

    //------------------------------------------------------------------------------------------------------------

    cout << endl << endl;

    double* v; // скорость звука
    v = new double[m];
    double avgv; // среднее значение скорости звука


    calculate_speed(v, u_m, m, f, d, a); // вычисление скорости с учетом дифр. поправки
    // d - длина образца
    // a - диаметр электрода
    // а - частота

    avgv = average(v, m - 1); // вычисялем среднее значение скорости звука
    cout << endl << "avgv=" << avgv << endl << endl;
    // посчитали и вывели скорость звука    

    //------------------------------------------------------------------------------------------------------------    

    double* alpha; // коэффициент затухания
    alpha = new double[m];
    double avgalpha; // среднее значение коэффициента затухания

    /*if (regime == "findmax") { // в режиме автокорреляции коэф. затухания не вычисляем
        calculate_attenuation(v, alpha, u_m, m, f, d, a);

        avgalpha = average(alpha, m - 1);

        cout << endl << "avgalpha=" << avgalpha << endl << endl;
        // посчитали и вывели коэффициент затухания
    }*/

    calculate_attenuation(v, alpha, u_m, m, f, d, a);

    avgalpha = average(alpha, m - 1);

    cout << endl << "avgalpha=" << avgalpha << endl << endl;
    // посчитали и вывели коэффициент затухания

        

    //------------------------------------------------------------------------------------------------------------


    file_out(u_m, v, alpha, m, d, f, a, filename, regime); //вывод в файл
    
    if (regime == "findmax") {
        firstmax[0] = u_m[0][0];
        firstmax[1] = u_m[1][0];
        cout << "firstmax=" << firstmax[0] << ", " << firstmax[1] << endl << endl;
    }

    delete u;
    delete u_max;
    delete u_m;
    delete v;
    delete alpha;


}
