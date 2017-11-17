#pragma once



namespace enhance
{
    template<typename DERIVED>
    struct Singleton
    {
        static DERIVED& getInstance()
        {
            if( INSTANCE == nullptr)
                INSTANCE = new DERIVED();
            return *INSTANCE;
        }
        
        static void destroyInstance()
        {
            delete INSTANCE;
            INSTANCE = nullptr;
        }
        
    protected:
        Singleton() { }
        
        virtual ~Singleton()
        {
            INSTANCE = nullptr;
        }
        
    private:
        static DERIVED* INSTANCE;
        
        Singleton(const Singleton&) = delete;
        Singleton& operator= (const Singleton&) = delete;
    };


    template<typename DERIVED> DERIVED* Singleton<DERIVED>::INSTANCE = nullptr;
}