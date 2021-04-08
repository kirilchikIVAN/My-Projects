#include <iostream>
#include <vector>
#include <memory>
#include <list>


template<size_t chunkSize>
class FixedAllocator {
private:
    size_t pool_size = 1 << 27;
    char* pool;
    size_t pointer;
public:

    FixedAllocator(): pool(new char[pool_size]), pointer(0) {}

    void* allocate(size_t count) {
        size_t delta = count * chunkSize;
        pointer += delta;
        return static_cast<void*>(pool + pointer - delta);
    }

    void* deallocate(void*) {
        return nullptr;
    }

    ~FixedAllocator() {
        delete [] pool;
    }
};


template<typename T>
class FastAllocator {
private:
    std::shared_ptr<FixedAllocator<1>> alloc_ptr;
    static const size_t MAX_BYTES = 128;

public:

    std::shared_ptr<FixedAllocator<1>> get_shared_ptr() const {
        return alloc_ptr;
    }

    typedef T value_type;
    FastAllocator(): alloc_ptr (std::make_shared<FixedAllocator<1>>()) {}

    template<typename U>
    FastAllocator<T>(const FastAllocator<U>& other): alloc_ptr(other.get_shared_ptr()) {}

    T* allocate(size_t count) {
        size_t size = count * sizeof(T);
        if (size <= MAX_BYTES) {
            return static_cast<T*>(alloc_ptr->allocate(size));
        }
        return static_cast<T*>(::operator new(size));
    }

    void* deallocate(T* ptr, size_t) {
        return alloc_ptr->deallocate(ptr);
    }

    template<typename... Args>
    void construct(T* ptr, const Args&... args) {
        new(ptr) T(args...);
    }

    void destroy(T* ptr) {
        ptr->~T();
    }

    ~FastAllocator<T>() = default;
};


template <typename T, typename U>
bool operator==(const FastAllocator<T>&, const FastAllocator<U>&) { return true; }

template <typename T, typename U>
bool operator!=(const FastAllocator<T>&, const FastAllocator<U>&) { return false; }




template<typename T, typename Allocator = std::allocator<T>>
class List {
public:

    struct Node {
        T value = T();
        Node* right = nullptr;
        Node* left = nullptr;

        Node() {}
        Node(const T& value): value(value) {}
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
        typename std::allocator_traits<Allocator>::template rebind_alloc<Node> NewAllocator;
        allocator = NewAllocator;
        fake = AllocTraits::allocate(allocator, 1);
        fake->right = fake;
        fake->left = fake;
    }

    List<T, Allocator>(size_t count, const T& value, const Allocator& alloc = Allocator()): T_allocator(alloc) {
        typename std::allocator_traits<Allocator>::template rebind_alloc<Node> NewAllocator;
        allocator = NewAllocator;
        fake = AllocTraits::allocate(allocator, 1);
        fake->left = fake;
        fake->right = fake;
        for (size_t i = 0; i < count; ++i) {
            push_back(value);
        }
    }

    List<T, Allocator>(size_t count, const Allocator& alloc = Allocator()): T_allocator(alloc) {
        typename std::allocator_traits<Allocator>::template rebind_alloc<Node> NewAllocator;
        allocator = NewAllocator;
        fake = AllocTraits::allocate(allocator, 1);
        fake->left = fake;
        fake->right = fake;
        len = count;
        for (size_t i = 0; i < count; ++i) {
            Node* new_node = AllocTraits::allocate(allocator, 1);
            AllocTraits::construct(allocator, new_node);
            new_node->left = fake->left;
            fake->left->right = new_node;
            new_node->right = fake;
            fake->left = new_node;
        }
    }

    Allocator get_allocator() const {
        return allocator;
    }

    List<T, Allocator>(const List<T, Allocator>& other) {
        allocator = AllocTraits::select_on_container_copy_construction(other.allocator);
        T_allocator = AllocTraits::select_on_container_copy_construction(other.T_allocator);
        fake = AllocTraits::allocate(allocator, 1);
        fake->left = fake;
        fake->right = fake;
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
        fake->left = new_node;
        new_node->right = fake;
        last_node->right = new_node;
        new_node->left = last_node;
    }

    void push_front(const T& value) {
        ++len;
        Node* new_node = AllocTraits::allocate(allocator, 1);
        AllocTraits::construct(allocator, new_node, value);
        Node* first_node = fake->right;
        fake->right = new_node;
        new_node->left = fake;
        first_node->left = new_node;
        new_node->right = first_node;
    }

    void pop_back() {
        if (empty()) {
            std::cerr << "Error: trying to pop_b empty list";
            return;
        }
        --len;
        Node* last_node = fake->left;
        Node* new_last_node = fake->left->left;
        fake->left = new_last_node;
        new_last_node->right = fake;

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
        fake->right = new_first_node;
        new_first_node->left = fake;

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

        left->right = new_node;
        new_node->left = left;
        right->left = new_node;
        new_node->right = right;
    }

    void erase(const_iterator it) {
        if (len == 0) {
            std::cerr << "Error: trying to erase in empty list/n";
            return;
        }
        --len;

        Node* left = it.ptr->left;
        Node* right = it.ptr->right;
        left->right = right;
        right->left = left;

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
