#ifndef ISIMULATIONS_H
#define ISIMULATIONS_H


class iSimulations
{
    public:
        iSimulations();
        virtual ~iSimulations();

        virtual bool one_freedom_degree()=0;
        virtual bool two_freedom_degree()=0;
        virtual bool tree_freedom_degree()=0;
        virtual bool save_data()=0;

    protected:

    private:
};

#endif // ISIMULATIONS_H
