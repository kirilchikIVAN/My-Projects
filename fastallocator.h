#include <iostream>
#include <vector>
#include <memory>
#include <list>
#include <queue>


template<size_t chunkSize>
class FixedAllocator {
private:
    static const size_t pool_size = 1 << 20;
    std::vector<char*> pool;
    char* pointer;
    std::queue<char*> free_ch;



    FixedAllocator() {
        pool.push_back(new char[pool_size]);
        free_ch.push(pool.back());
        pointer = pool.back();
    }

public:

    void* allocate() {
        char* chunk = free_ch.front();
        if (chunk == pointer) {
            pointer += chunkSize;
            if (pointer >= pool.back() + pool_size) {
                pool.push_back(new char[pool_size]);
                pointer = pool.back();
            }
            free_ch.push(pointer);

        }
        free_ch.pop();
        return static_cast<void*>(chunk);
    }

    void deallocate(void* ptr) {
        free_ch.push(static_cast<char*>(ptr));
    }

    static FixedAllocator& getInstance() {
        static FixedAllocator* self;
        if (self == nullptr) self = new FixedAllocator();
        return *self;
    }

    ~FixedAllocator() {
        for (auto* p : pool) {
            delete [] p;
        }
    }
};


template<typename T>
class FastAllocator {
private:

public:


    using value_type = T;

    FastAllocator() = default;

    template<typename U>
    FastAllocator<value_type>(const FastAllocator<U>&) {}

    FastAllocator& operator=(const FastAllocator&) = default;

    value_type* allocate(size_t count) {
        size_t size = count * sizeof(T);
        void* ptr;

        if (size <= 16) ptr = FixedAllocator<16>::getInstance().allocate();
        else if (size <= 24) ptr = FixedAllocator<24>::getInstance().allocate();
        else if (size <= 32) ptr = FixedAllocator<32>::getInstance().allocate();
        else if (size <= 64) ptr = FixedAllocator<64>::getInstance().allocate();
        else ptr = ::operator new(size);
        return static_cast<T*>(ptr);
    }

    void deallocate(value_type* ptr, size_t count) {
        size_t size = sizeof(T) * count;

        if (size <= 16) FixedAllocator<16>::getInstance().deallocate(ptr);
        else if (size <= 24) FixedAllocator<24>::getInstance().deallocate(ptr);
        else if (size <= 32) FixedAllocator<32>::getInstance().deallocate(ptr);
        else if (size <= 64) FixedAllocator<64>::getInstance().deallocate(ptr);
        else ::operator delete(ptr, size);
    }

    ~FastAllocator<T>() = default;
};

template <typename T, typename U>
bool operator==(const FastAllocator<T>&, const FastAllocator<U>&) { return true; }

template <typename T, typename U>
bool operator!=(const FastAllocator<T>& first, const FastAllocator<U>& second) { return !(first == second); }




template<typename T, typename Allocator = std::allocator<T>>
class List {
public:

    struct Node {
        T value = T();
        Node* right = nullptr;
        Node* left = nullptr;

        Node() {}
        Node(const T& value): value(value) {}

        void link(Node* other) {
            right = other;
            other->left = this;
        }
    };

    typename std::allocator_traits<Allocator>::template rebind_alloc<Node> allocator;
    Allocator T_allocator;
    using AllocTraits = std::allocator_traits<typename std::allocator_traits<Allocator>::template rebind_alloc<Node>>;

    Node* fake = nullptr;
    size_t len = 0;

    bool empty() const {
        return len == 0;
    }

    void swap(List<T, Allocator>& other) {
        std::swap(fake, other.fake);
        std::swap(len, other.len);
    }

public:

    void show() const {
        Node* tmp_node = fake->right;
        while (tmp_node != fake) {
            tmp_node = tmp_node->right;
        }
    }

    explicit List<T, Allocator>(const Allocator& alloc = Allocator()): T_allocator(alloc), len(0)  {
        typename std::allocator_traits<Allocator>::template rebind_alloc<Node> FixedAllocator;
        allocator = FixedAllocator;
        fake = AllocTraits::allocate(allocator, 1);
        fake->link(fake);
    }

    List<T, Allocator>(size_t count, const T& value, const Allocator& alloc = Allocator()): T_allocator(alloc) {
        typename std::allocator_traits<Allocator>::template rebind_alloc<Node> FixedAllocator;
        allocator = FixedAllocator;
        fake = AllocTraits::allocate(allocator, 1);
        fake->link(fake);
        for (size_t i = 0; i < count; ++i) {
            push_back(value);
        }
    }

    List<T, Allocator>(size_t count, const Allocator& alloc = Allocator()): T_allocator(alloc) {
        typename std::allocator_traits<Allocator>::template rebind_alloc<Node> FixedAllocator;
        allocator = FixedAllocator;
        fake = AllocTraits::allocate(allocator, 1);
        fake->link(fake);
        len = count;
        for (size_t i = 0; i < count; ++i) {
            Node* new_node = AllocTraits::allocate(allocator, 1);
            AllocTraits::construct(allocator, new_node);
            fake->left->link(new_node);
            new_node->link(fake);
        }
    }

    Allocator get_allocator() const {
        return allocator;
    }

    List<T, Allocator>(const List<T, Allocator>& other) {
        allocator = AllocTraits::select_on_container_copy_construction(other.allocator);
        T_allocator = AllocTraits::select_on_container_copy_construction(other.T_allocator);
        fake = AllocTraits::allocate(allocator, 1);
        fake->link(fake);
        auto it = other.begin();
        while (it != other.fake) {
            insert(fake, *it);
            ++len;
            ++it;
        }
        len = other.len;
    }

    List& operator=(const List& other) {
        if (this == &other) return *this;

        if (AllocTraits::propagate_on_container_copy_assignment::value && allocator != other.allocator) {
            allocator = other.allocator;
            T_allocator = other.T_allocator;
        }
        while (!empty()) {
            pop_back();
        }
        List<T, Allocator> copy(other);
        swap(copy);
        return *this;
    }

    size_t size() const {
        return len;
    }

    void push_back(const T& value) {
        ++len;
        Node* new_node = AllocTraits::allocate(allocator, 1);
        AllocTraits::construct(allocator, new_node, value);
        Node* last_node = fake->left;
        new_node->link(fake);
        last_node->link(new_node);
    }

    void push_front(const T& value) {
        ++len;
        Node* new_node = AllocTraits::allocate(allocator, 1);
        AllocTraits::construct(allocator, new_node, value);
        Node* first_node = fake->right;
        fake->link(new_node);
        new_node->link(first_node);
    }

    void pop_back() {
        if (empty()) {
            std::cerr << "Error: trying to pop_b empty list";
            return;
        }
        --len;
        Node* last_node = fake->left;
        Node* new_last_node = fake->left->left;
        new_last_node->link(fake);

        AllocTraits::destroy(allocator, last_node);
        AllocTraits::deallocate(allocator, last_node, 1);

    }

    void pop_front() {
        if (empty()) {
            std::cerr << "Error: trying to pop_b empty list";
            return;
        }
        --len;
        Node* first_node = fake->right;
        Node* new_first_node = fake->right->right;
        fake->link(new_first_node);

        AllocTraits::destroy(allocator, first_node);
        AllocTraits::deallocate(allocator, first_node, 1);
    }

    template<bool is_const>
    class iterator_base {
    private:
        Node* ptr;
        iterator_base<is_const>(Node* other): ptr(other) {}
    public:
        friend class List;

        using pointer = typename std::conditional_t<is_const, const T*, T*>;
        using reference = typename std::conditional_t<is_const, const T&, T&>;
        using value_type = typename std::conditional_t<is_const, const T, T>;

        using difference_type = std::ptrdiff_t;
        using iterator_category = std::bidirectional_iterator_tag;

        iterator_base& operator++() {
            ptr = ptr->right;
            return *this;
        }

        iterator_base& operator--() {
            ptr = ptr->left;
            return *this;
        }

        iterator_base operator++(int) {
            iterator_base copy = *this;
            ptr = ptr->right;
            return copy;
        }

        iterator_base operator--(int) {
            iterator_base copy = *this;
            ptr = ptr->left;
            return copy;
        }

        bool operator==(const iterator_base& other) const {
            return ptr == other.ptr;
        }

        bool operator!=(const iterator_base& other) const {
            return !(operator==(other.ptr));
        }

        reference operator*() const {
            return ptr->value;
        }

        pointer operator->() const {
            return &ptr->value;
        }

        operator iterator_base<true>() const {
            return iterator_base<true>(ptr);
        }
    };

    using iterator = iterator_base<false>;
    using const_iterator = iterator_base<true>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    iterator begin() { return iterator(fake->right); }

    const_iterator begin() const { return const_iterator(fake->right); }

    iterator end() { return iterator(fake); }

    const_iterator end() const { return const_iterator(fake); }

    const_iterator cbegin() const { return const_iterator(fake->right); }

    const_iterator cend() const { return const_iterator(fake); }

    reverse_iterator rbegin() { return reverse_iterator(fake); }

    const_reverse_iterator rbegin() const { return const_reverse_iterator(fake); }

    reverse_iterator rend() { return reverse_iterator(fake->right); }

    const_reverse_iterator rend() const { return const_reverse_iterator(fake->right); }

    const_reverse_iterator crbegin() { return const_reverse_iterator(fake); }

    const_reverse_iterator crend() { return const_reverse_iterator(fake->right); }

    void insert(const_iterator it, const T& value) {
        ++len;
        Node* new_node = AllocTraits::allocate(allocator, 1);
        AllocTraits::construct(allocator, new_node, value);
        Node* left = it.ptr->left;
        Node* right = it.ptr;

        left->link(new_node);
        new_node->link(right);
    }

    void erase(const_iterator it) {
        if (len == 0) {
            std::cerr << "Error: trying to erase in empty list/n";
            return;
        }
        --len;

        Node* left = it.ptr->left;
        Node* right = it.ptr->right;
        left->link(right);

        AllocTraits::destroy(allocator, it.ptr);
        AllocTraits::deallocate(allocator, it.ptr, 1);
    }

    ~List() {
        if (!empty()) {
            iterator it = begin();
            while (it.ptr != fake) {
                iterator to_destroy = it++;
                AllocTraits::destroy(allocator, to_destroy.ptr);
                AllocTraits::deallocate(allocator, to_destroy.ptr, 1);
            }
        }
        AllocTraits::deallocate(allocator, fake, 1);
    }
};
