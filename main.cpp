#include <conio.h>
#include <windows.h>

#include "Display.h" // namespace
#include "iSimulations.h" // interface
#include "Simulations_v1.h" // singleton

#define ESC 27

using namespace Display;

int main()
{
    Simulations_v1& singleSimulation=Simulations_v1::get_singleSimulation();
    char selection;
    do
    {
        system("cls");
        menu_display();
        selection=getch();

        if(selection=='1')
        {
            system("cls");
            simulation_selection_display();
            selection=getch();
            if(selection=='1')
            {
                system("cls");
                singleSimulation.one_freedom_degree();
            }
            else if(selection=='2')
            {
                system("cls");
                singleSimulation.two_freedom_degree();
            }
            else if(selection=='3')
            {
                system("cls");
                singleSimulation.tree_freedom_degree();
            }
        }
        if(selection=='2')
        {
            singleSimulation.save_data();
            system("pause");
        }

    }while(selection!=ESC);

    return 0;
}
