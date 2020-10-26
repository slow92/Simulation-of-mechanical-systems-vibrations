#include "Display.h"

#include<iostream>

using std::cout;
using std::endl;

namespace Display
{
    void menu_display()
    {
        cout<<"[1] New simulation"<<endl;
        cout<<"[2] Save simulation results"<<endl;
        cout<<"ESC to exit"<<endl;
    }

    void simulation_selection_display()
    {
        cout<<"[1] One_freedom_degree"<<endl;
        cout<<"[2] Two_freedom_degree"<<endl;
        cout<<"[3] Tree_freedom_degree"<<endl;
    }
}
