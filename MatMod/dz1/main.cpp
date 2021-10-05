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
int main() {
    std::string line;
    int xn,yn,fx,fy;
    int finde_value=1;
    double min_d=1;
    bool isFirstLine= true;
    std::ifstream in("H:\\Ucheba\\Poly_3_kurs\\3_Kurs_MexMatMod\\MatMod\\dz1\\files_for_dz\\in.txt"); // окрываем файл для чтения
    if (in.is_open())
    {
        while (getline(in, line))
        {
            string *values=StringToMass(line,' ',2);
            int x=atoi( values[0].c_str() );
            int y=atoi( values[1].c_str() );
            if(isFirstLine){
                isFirstLine= false;
                xn=(-1)*x;
                yn=(-1)*y;
            }
            else{
                //<0 - справа
                if(((xn-0)*(y-0)-(yn-0)*(x-0))<=0){
                    double _cos=(xn*x+yn*y)/(sqrt(pow(x,2)+pow(y,2))*sqrt(pow(xn,2)+pow(yn,2)));
                    std::cout << _cos << std::endl;
                    if(finde_value-_cos<min_d){
                        min_d=finde_value-_cos;
                        fx=x;
                        fy=y;
                    }
                }
            }
            std::cout << line << std::endl;
            delete[] values;
        }
    }
    in.close();     // закрываем файл
    std::cout << "Antworten:" << std::endl;
    std::cout << "delta="+to_string(min_d)+" x="+to_string(fx)+" y="+to_string(fy) << std::endl;
    return 0;
}
