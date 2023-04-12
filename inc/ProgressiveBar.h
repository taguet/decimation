#ifndef PROGRESSIVEBAR_H
#define PROGRESSIVEBAR_H

#include <ostream>
 
class ProgressiveBar
{
    public:
        ProgressiveBar(std::ostream &os, int const& length=30);
        void init();
        void update(int const& i, int const& n);
         
    private:
        int m_length;
        int m_mark;
        int m_progress;
        std::ostream &m_os;
};
 
inline void ProgressiveBar::update(int const& i, int const& n)
{
    int p = int(((float(i/(n-1)))*m_length));
    m_progress=p;
    if (m_progress>m_mark)
    {
        m_mark=m_progress;
        m_os << ".";
    }
    m_os.flush();
}
 
#endif//PROGRESSIVEBAR_H
