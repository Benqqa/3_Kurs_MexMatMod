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
int main() {
    std::string line;
    int xn,yn,fx_r,fy_r,fx_l,fy_l;
    int finde_value=1;
    double min_d_r=1;
    double min_d_l=1;
    bool isFirstLine= true;
    bool isTwiceLine_r= true,isTwiceLine_l= true;
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

                //растояние
                double dist;
                if(xn !=0 && yn !=0){
                    dist=(abs(-1*(x)/xn+y/yn)/sqrt(pow(1/xn,2)+pow(1/yn,2)));
                   // dist=(abs((yn-0)*x-(xn-0)*y+xn*0-yn*0))/(sqrt(pow(xn-0,2)+pow(yn-0,2))); //вектор как прямая заданная 2мя точками
                }
                else{
                    if(xn ==0){
                        dist= x;
                        //dist=sqrt(pow(x-xn,2)+pow(y-yn,2)); //как растояние между 2мя точками
                    }
                    else{
                        dist=y;
                    }
                }
                //<0 - справа"Left:
                if(((xn-0)*(y-0)-(yn-0)*(x-0))<=0){
                    std::cout << "Right, Dist="+to_string(dist) << std::endl;
                    if(isTwiceLine_r){
                        isTwiceLine_r= false;
                        min_d_r=dist;
                    }
                    else{
                        if(dist<min_d_r){
                            min_d_r=dist;
                            fx_r=x;
                            fy_r=y;
                        }
                    }
                }//>0 - слева
                else{
                    std::cout << "Left, Dist="+to_string(dist) << std::endl;
                    if(isTwiceLine_l){
                        isTwiceLine_l= false;
                        min_d_l=dist;
                    }
                    else{
                        if(dist<min_d_l){
                            min_d_l=dist;
                            fx_l=x;
                            fy_l=y;
                        }
                    }
                }
                //
                //угол
                /*
                double _cos=(xn*x+yn*y)/(sqrt(pow(x,2)+pow(y,2))*sqrt(pow(xn,2)+pow(yn,2)));
                std::cout << _cos << std::endl;
                if(finde_value-_cos<min_d){
                    min_d=finde_value-_cos;
                    fx=x;
                    fy=y;
                }
                 */
                //
            }
            std::cout << line << std::endl;
            delete[] values;
        }
    }
    in.close();     // закрываем файл
    std::cout << "Antworten:" << std::endl;
    std::cout << "Left: delta="+to_string(min_d_l)+" x="+to_string(fx_l)+" y="+to_string(fy_l) << std::endl;
    std::cout << "Right: delta="+to_string(min_d_r)+" x="+to_string(fx_r)+" y="+to_string(fy_r) << std::endl;

    return 0;
}
