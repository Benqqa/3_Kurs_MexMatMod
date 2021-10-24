#include <iostream>
#include <fstream>
#include <iomanip>
#include <valarray>

using namespace std;
string * StringToMass(string base_str, char delim, int size){

    size_t pos = 0;
    size_t base_str_size = base_str.size();
    //size_t delim_size = delim.size();
    size_t delim_size=1;
    string temp;
    int i=0;
    string *mass_str = new string[size];
    while (pos < base_str_size) {
        temp = temp.assign(base_str, pos, base_str.find(delim, pos) - pos);
        if (temp.size() > 0)  // проверка на пустую строку при необходимости
        {
            // cout << temp << endl;
            mass_str[i] = temp;
            i++;
            //cout << i << endl;
        }
        pos += temp.size() + delim_size;

    }

    return mass_str;
}
std::tuple<int, int> step(int dividend, int divisor) {
    return  std::make_tuple(dividend / divisor, dividend % divisor);
}
void polet(int index, int mass_size,string *mass_stolbov,int napr, double g, double y_0,double x_0, double x_d, double v_x, double v_y, double v_0, int nomer_stolba){
    std::cout << "index = "+to_string(index) << std::endl;
    if(index == 0){// первый вход - наполним массив столбов до первого пападания

        std::string line;
        std::ifstream in("H:\\Ucheba\\Poly_3_kurs\\3_Kurs_MexMatMod\\MatMod\\dz2\\files_for_dz\\in.txt"); // окрываем файл для чтения
        std::ifstream in2("H:\\Ucheba\\Poly_3_kurs\\3_Kurs_MexMatMod\\MatMod\\dz2\\files_for_dz\\in.txt"); // окрываем файл для чтения
        //std::cout << line << std::endl;
        if (in.is_open())
        {

            int m=0;
            while (getline(in, line)) // перебераем столбы
            {
                if(line.length() != 0){
                    m++;
                }
            }

            string *mass_stolbov = new string[m];

            int k= 0;
            while (getline(in2, line)) // перебераем столбы
            {

            //кусок кода надо формить красиво
                k++;

                if(k == 1){
                    mass_stolbov[k-1]="0 "+line;
                    string *values_start=StringToMass("0 "+line,' ',2);
                    y_0=atof( values_start[1].c_str() );
                    x_0=atof( values_start[0].c_str() );
                    std::cout << values_start[0]+" "+ values_start[1] << std::endl;
                    delete[] values_start;
                    std::cout << "skip"+line << std::endl;
                    continue;
                }
                if(k == 2){
                    string *values_start=StringToMass(line,' ',2);
                    x_d=0; // коэф смещения
                    v_x=atof( values_start[0].c_str() );
                    v_y=atof( values_start[1].c_str() );
                    v_0= sqrt(pow(v_x,2)+pow(v_y,2));
                    delete[] values_start;
                    std::cout << "skip"+line << std::endl;
                    continue;
                }
                std::cout << line << std::endl;
                std::cout << to_string(v_x)+" "+ to_string(v_y)+" "+ to_string(v_0) << std::endl;
                std::cout << to_string(x_0)+" "+ to_string(y_0) << std::endl;
                std::cout << "null "+mass_stolbov[0] << std::endl;
                string *values=StringToMass(line,' ',2);
                double x=atof( values[0].c_str() );
                double y=atof( values[1].c_str() );
                std::cout << to_string(x)+" "+ to_string(y) << std::endl;
                double new_y;
                /*if( napr==1 && x > x_0){ //отбор столбиков от направления // летит в право // а нафига это вообще?????
                    napr = -1;*/
                    //проверка столбика
                    new_y=y_0+napr*(v_y/v_x)*(x-x_0-2*x_d)-g*(pow(x-x_0-2*x_d,2)/(2*pow(v_0,2)));
                    std::cout << "Heigth = "+to_string(new_y) << std::endl;
                    if(x==x_0){
                        std::cout <<"VSE mi teper letim nazad!!! uuuhhuuuu x_0=x" << std::endl;
                        std::cout <<"ZONA = 0" << std::endl;
                        break;
                        break;
                    }
                    if(new_y <= 0){
                        std::cout << "Upalo v 0" << std::endl;
                        std::cout <<"ZONA = "+to_string(nomer_stolba) << std::endl;
                        break;
                    }
                    if(new_y<=y){ // попал

                        mass_stolbov[k-2]=line;
                        std::cout << "POPALO!!!! v stolb x= "+to_string(x)+" y= "+to_string(y) << std::endl;
                        x_d=x-x_d; //смещение
                        //рекурсим
                        index++;
                        //nomer_stolba=nomer_stolba+napr*k-1;
                        napr=-1*napr;
                        //------
                        //тут надо развернуть массив координат столбов
                        //-------
                        std::cout << k << std::endl;
                        if(k-1!= 2){

                            string *mass_stolbov_2 = new string[k-1];
                            std::cout << "mass_stolbov_2" << std::endl;
                            for(int i=0; i<k-1; i++){
                                mass_stolbov_2[i]=mass_stolbov[i];
                                std::cout <<mass_stolbov_2[i] << std::endl;
                            }
                            delete[] mass_stolbov;
                            std::cout << mass_stolbov_2->length() << std::endl;
                            string *mass_stolbov_2_1 = new string[k-1];
                            std::cout <<"perevorot" << std::endl;
                            for(int i=0; i<k-1; i++){
                                std::cout <<mass_stolbov_2[k-2-i] << std::endl;
                                mass_stolbov_2_1[i]=mass_stolbov_2[k-2-i];

                            }
                            for(int i=0; i<k-1; i++){
                                mass_stolbov_2[i]=mass_stolbov_2_1[i];
                            }
                            delete[] mass_stolbov_2_1;
                            mass_size=k-1;
                            polet(index,mass_size,mass_stolbov_2,napr, g, y_0, x_0,  x_d,  v_x,  v_y,  v_0, nomer_stolba);
                        }
                        else{
                            std::cout <<"VSE!!! megdu 2 stolbami" << std::endl;
                            std::cout <<"ZONA = "+to_string(nomer_stolba) << std::endl;
                            break;
                        }
                    }
                    else{ // получается не попал)))

                        std::cout << "NE POPALO(((( x= "+to_string(x)+" y= "+to_string(y) << std::endl;
                        nomer_stolba+=napr*1;
                        //------------
                        //тут надо слхранить координаты в массив
                        mass_stolbov[k-2]=line;
                        //--------------
                        if(k-1 == m-1){
                            std::cout << "Vse pezda! - pereletelo" << std::endl;
                            std::cout <<"ZONA = "+to_string(nomer_stolba) << std::endl;
                        }
                    }
                delete[] values;
            }
        }
    }
//-----------------------------------------------------------------------------------------------
    else{
    //кусок кода надо формить красиво
        int k=0;
        string line;
        double new_y;
       // std::cout << mass_stolbov[0] << std::endl;
        std::cout << "mass_stolbov/*" << std::endl;
        std::cout << mass_size << std::endl;
       for(int i=0; i<mass_size; i++){
           std::cout << mass_stolbov[i] << std::endl;
       }
        std::cout << "*/" << std::endl;
        string *mass_stolbov_1 = new string[mass_size];
        while (true){
            //mass_stolbov[k] доставать иксы и делать цирк

            line =mass_stolbov[k];
            k++;
            if(k == 1){
                mass_stolbov_1[k-1]=line;
                continue;
            }
            std::cout << line << std::endl;
            string *values=StringToMass(line,' ',2);
            double x=atof( values[0].c_str() );
            double y=atof( values[1].c_str() );
            /*if( napr==1 && x > x_0){ //отбор столбиков от направления // летит в право // а нафига это вообще?????
                napr = -1;*/
            //проверка столбика
            new_y=y_0+napr*(v_y/v_x)*(x-x_0-2*x_d)-g*(pow(x-x_0-2*x_d,2)/(2*pow(v_0,2)));
            std::cout << "Heigth = "+to_string(new_y) << std::endl;
            if(x==x_0){
                std::cout <<"VSE mi teper letim nazad!!! uuuhhuuuu x_0=x" << std::endl;
                std::cout <<"ZONA = 0" << std::endl;
                break;
                break;
            }
            if(new_y <= 0){
                std::cout << "Upalo v 0" << std::endl;
                std::cout <<"ZONA = "+to_string(nomer_stolba) << std::endl;
                break;
            }
            if(new_y<=y){ // попал
                mass_stolbov_1[k-1]=line;
                std::cout << "POPALO!!! x= "+to_string(x)+" y= "+to_string(y) << std::endl;
                x_d=x-x_d; //смещение
                //рекурсим
                index++;
                //nomer_stolba=nomer_stolba+napr*k;
                std::cout <<"ZONA = "+to_string(nomer_stolba) << std::endl;
                napr=-1*napr;
                //------
                //тут разворот массива координат столбов
                std::cout << k << std::endl;
                if(k!=2){

                    string *mass_stolbov_2 = new string[k];
                    for(int i=0; i<k; i++){
                        mass_stolbov_2[i]=mass_stolbov_1[i];
                        std::cout <<mass_stolbov_2[i] << std::endl;
                    }
                    delete[] mass_stolbov;
                    delete[] mass_stolbov_1;
                    std::cout << "len: "+to_string(mass_stolbov_2->length()) << std::endl;
                    string *mass_stolbov_2_1 = new string[k];
                    std::cout <<"perevorot" << std::endl;
                    for(int i=0; i<k; i++){
                        std::cout <<mass_stolbov_2[k-1-i] << std::endl;
                        mass_stolbov_2_1[i]=mass_stolbov_2[k-1-i];

                    }
                    for(int i=0; i<k; i++){
                        mass_stolbov_2[i]=mass_stolbov_2_1[i];
                    }
                    delete[] mass_stolbov_2_1;
                    mass_size=k;
                    //-------
                    polet(index,mass_size,mass_stolbov_2,napr, g, y_0, x_0,  x_d,  v_x,  v_y,  v_0, nomer_stolba);
                }
                else{
                    std::cout <<"VSE!!! megdu 2 stolbami" << std::endl;
                    std::cout <<"ZONA = "+to_string(nomer_stolba) << std::endl;
                    break;
                }
            }
            else{ // получается не попал)))
                //------------
                //слхраним координаты в массив
                if(x!=x_0){
                    std::cout << "NE POPAPO(((( x= "+to_string(x)+" y= "+to_string(y) << std::endl;
                    nomer_stolba+=napr*1;
                    mass_stolbov_1[k-1]=line;
                    //--------------
                }
                else{

                    std::cout <<"VSE mi teper letim nazad!!! uuuhhuuuu x_0=x" << std::endl;
                    std::cout <<"ZONA = 0" << std::endl;
                    break;
                    break;
                }


            }
            delete[] values;
        }
    }
}
int main(int argc, char** argv) {
    if(argc == 2){
        // есть один агрумент
        // в argv[1] содержится строка с первым агрументом (имя файла)
        std::cout << "1st argument: "<< argv[1] << std::endl;
    }else{
        // аргументов нет или их больше чем мы ожидаем
    }

    //на вход данные о полете
    double g=9.81;
    double y_0=1;
    double x_0=0;
    double x_d=0; // коэф смещения
    double v_x=3;
    double v_y=1;
    double v_0= sqrt(pow(v_x,2)+pow(v_y,2));
    int napr = 1; //напрво - true; налево - false
    int index = 0;
    int mass_size=0;
    int nomer_stolba=0;
    std::string line;
    std::ifstream in("H:\\Ucheba\\Poly_3_kurs\\3_Kurs_MexMatMod\\MatMod\\dz2\\files_for_dz\\in.txt"); // окрываем файл для чтения

    if (in.is_open()) {
        int m = 0;
        while (getline(in, line)) // перебераем столбы
        {
            std::cout << line << std::endl;
            if (line.length() != 0) {
                m++;
            }
        }
        mass_size=m;
        string *mass_stolbov = new string[m];
        polet(index,mass_size, mass_stolbov, napr, g, y_0, x_0, x_d, v_x, v_y, v_0,nomer_stolba);
    }
    //
    return 0;
}
