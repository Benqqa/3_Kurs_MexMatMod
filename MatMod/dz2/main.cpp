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
void polet(int index,string *mass_stolbov,int napr, int g, double y_0,double x_0, double x_d, double v_x, double v_y, double v_0, int nomer_stolba){
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

            string *mass_stolbov = new string[m-1];

            int k= 0;
            while (getline(in2, line)) // перебераем столбы
            {

            //кусок кода надо формить красиво
                k++;
                if(k == 1){
                    continue;
                }
                std::cout << line << std::endl;
                string *values=StringToMass(line,' ',2);
                double x=atof( values[0].c_str() );
                double y=atof( values[1].c_str() );
                double new_y;
                /*if( napr==1 && x > x_0){ //отбор столбиков от направления // летит в право // а нафига это вообще?????
                    napr = -1;*/
                    //проверка столбика
                    new_y=y_0+napr*(v_y/v_x)*(x-x_0-2*x_d)-g*(pow(x-x_0-2*x_d,2)/(2*pow(v_0,2)));
                    std::cout << "Heigth = "+to_string(new_y) << std::endl;
                    if(new_y == 0){
                        std::cout << "Upalo v 0" << std::endl;
                        break;
                    }
                    if(new_y<=y){ // попал
                        mass_stolbov[k-2]=line;
                        std::cout << "POPALO!!!! v stolb x= "+to_string(x)+" y= "+to_string(y) << std::endl;
                        x_d=x-x_d; //смещение
                        //рекурсим
                        index++;
                        nomer_stolba=nomer_stolba+napr*k-1;
                        napr=-1*napr;
                        //------
                        //тут надо развернуть массив координат столбов
                        //-------
                        std::cout << k << std::endl;
                        if(k!= 2){

                            string *mass_stolbov_2 = new string[k-1];
                            for(int i=0; i<k-1; i++){
                                mass_stolbov_2[i]=mass_stolbov[i];
                                std::cout <<mass_stolbov_2[i] << std::endl;
                            }
                            delete[] mass_stolbov;
                            std::cout << mass_stolbov_2->length() << std::endl;
                            string *mass_stolbov_2_1 = new string[k-1];
                            std::cout <<"perevorot" << std::endl;
                            for(int i=0; i<k-1; i++){
                                //std::cout <<mass_stolbov_2[i] << std::endl;
                                //std::cout <<"i= "+to_string(k-3-i) << std::endl;
                                std::cout <<mass_stolbov_2[k-2-i] << std::endl;
                                mass_stolbov_2_1[i]=mass_stolbov_2[k-2-i];

                            }
                            for(int i=0; i<k-1; i++){
                                // std::cout <<mass_stolbov_2[i] << std::endl;
                                mass_stolbov_2[i]=mass_stolbov_2_1[i];
                            }
                            delete[] mass_stolbov_2_1;
                            polet(index,mass_stolbov_2,napr, g, y_0, x_0,  x_d,  v_x,  v_y,  v_0, nomer_stolba);
                        }
                        else{
                            std::cout <<"VSE!!! megdu 2 stolbami" << std::endl;
                            std::cout <<"ZONA = "+to_string(nomer_stolba) << std::endl;
                            break;
                        }
                    }
                    else{ // получается не попал)))
                        std::cout << "NE POPALO(((( x= "+to_string(x)+" y= "+to_string(y) << std::endl;
                        //------------
                        //тут надо слхранить координаты в массив
                        mass_stolbov[k-2]=line;
                        //--------------

                    }
                delete[] values;
            }
        }
    }
    else{
    //кусок кода надо формить красиво
        int k=0;
        string line;
        double new_y;
       // std::cout << mass_stolbov[0] << std::endl;
        while (true){
            //mass_stolbov[k] доставать иксы и делать цирк
            line =mass_stolbov[k];
            k++;
            if(k == 1){
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
            if(new_y == 0){
                std::cout << "Upalo v 0" << std::endl;
                break;
            }
            if(new_y<=y){ // попал
                mass_stolbov[k-2]=line;
                std::cout << "POPALO!!! x= "+to_string(x)+" y= "+to_string(y) << std::endl;
                x_d=x-x_d; //смещение
                //рекурсим
                index++;
                nomer_stolba=nomer_stolba+napr*k-1;
                napr=-1*napr;
                //------
                //тут надо развернуть массив координат столбов
                std::cout << k << std::endl;
                if(k!=2){

                    string *mass_stolbov_2 = new string[k-1];
                    for(int i=0; i<k-1; i++){
                        mass_stolbov_2[i]=mass_stolbov[i];
                        std::cout <<mass_stolbov_2[i] << std::endl;
                    }
                    delete[] mass_stolbov;
                    std::cout << "len: "+to_string(mass_stolbov_2->length()) << std::endl;
                    string *mass_stolbov_2_1 = new string[k-1];
                    std::cout <<"perevorot" << std::endl;
                    for(int i=0; i<k-1; i++){
                        //std::cout <<mass_stolbov_2[i] << std::endl;
                        //std::cout <<"i= "+to_string(k-3-i) << std::endl;
                        std::cout <<mass_stolbov_2[k-3-i] << std::endl;
                        mass_stolbov_2_1[i]=mass_stolbov_2[k-2-i];

                    }
                    for(int i=0; i<k-1; i++){
                        // std::cout <<mass_stolbov_2[i] << std::endl;
                        mass_stolbov_2[i]=mass_stolbov_2_1[i];
                    }
                    delete[] mass_stolbov_2_1;
                    //-------
                    polet(index,mass_stolbov_2,napr, g, y_0, x_0,  x_d,  v_x,  v_y,  v_0, nomer_stolba);
                }
                else{
                    std::cout <<"VSE!!! megdu 2 stolbami" << std::endl;
                    std::cout <<"ZONA = "+to_string(nomer_stolba) << std::endl;
                    break;
                }
            }
            else{ // получается не попал)))
                //------------
                //тут надо слхранить координаты в массив
                std::cout << "NE POPAPO(((( x= "+to_string(x)+" y= "+to_string(y) << std::endl;
                mass_stolbov[k-2]=line;
                //--------------

            }
            delete[] values;
        }
    }
}
int main() {

    //на вход данные о полете
    int g=10;
    double y_0=1;
    double x_0=0;
    double x_d=0; // коэф смещения
    double v_x=3;
    double v_y=1;
    double v_0= sqrt(pow(v_x,2)+pow(v_y,2));
    v_0=3.1;
    int napr = 1; //напрво - true; налево - false
    int index = 0;
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

        string *mass_stolbov = new string[m-1];
        polet(index, mass_stolbov, napr, g, y_0, x_0, x_d, v_x, v_y, v_0,nomer_stolba);
    }
    //
    /*
    //данные о столбиках
    double x = 1;
    double y = 2;
    //отбор столбиков от направления
    if( napr==1 && x > x_0){

    }
    if( napr != 1 && x < x_0){

    }

    //
    //проверка столбика
    double new_y=y_0+napr*(v_y/v_x)*(x-x_0-2*x_d)-g*(pow(x-x_0-2*x_d,2)/(2*pow(v_0,2)));
    if(new_y<=y){
        x_d=x-x_d; //смещение
        //рекурсим

    }
    std::string line;
    std::ifstream in("H:\\Ucheba\\Poly_3_kurs\\3_Kurs_MexMatMod\\MatMod\\dz2\\files_for_dz\\in.txt"); // окрываем файл для чтения
    if (in.is_open())
    {
        while (getline(in, line))
        {
            std::cout << line << std::endl;

        }
    }*/
    return 0;
}
