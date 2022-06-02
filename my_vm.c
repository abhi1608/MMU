#include "my_vm.h"

//Helper Functions
void split_virtual_bits();

void intialize_PT();

static void set_bit(char *bit_map, int index);

static int get_bit(char *bit_map, int index);

static void reset_bit(char *bit_map, int index);

bool init = false;


#define page_array_size (PGSIZE/sizeof(unsigned long))

// Number of physical pages
unsigned long long physical_page_num = (MEMSIZE)/(PGSIZE);

//Number of virtual pages
unsigned long long virtual_page_num = (MAX_MEMSIZE)/(PGSIZE);


// To have a physical abstraction to pages, we divide into an array of unsigned longs.
// These array  had number of longs that fit into a single page.
struct page {
    unsigned long arr[page_array_size];
};

// Pointer to the actual physical memory, since 
// NOTE : memory is accessed as pages.
struct page *physical_memory;

// Physical bit map
char *p_bit_map;

// Virtual bit map
char *v_bit_map;

//pde_t * ptrToPD = NULL;

// Globals used to store various bits that are extracted from Virtual Address
int inner_level_bits = 0;
int outer_level_bits = 0;
int offset_bits = 0;

// To store tlb misses
int num_tlb_misses = 0;

//To store tlb lookups
int num_tlb_lookups = 0;

// Mutex for CS
pthread_mutex_t mutex;


// This function is responsible for splitting the 32 bit virtual address.
// It calculates the offset bits, bits needed for inner page table and
// number of bits needed for outer bits.
void split_virtual_bits()
{	
    //Get the bits for page offset
    offset_bits = log2(PGSIZE);

    //Get the bits of VPN
    int vpn_bits = (32 - offset_bits);

    // Calculating the inner level bits by considering what is taught in class
    // that page_size = (no_of_entries) * sizeof_pte
    inner_level_bits = log2((PGSIZE)/sizeof(pte_t));

    // outer level bits are the remaining bits 
    outer_level_bits = (vpn_bits - inner_level_bits);


}

void intialize_PT()
{
    // Reserve last page of physical memory for the page directory. 
    // This will have 2 ^ outer_level_bits number of 2nd level page tables 
	int index = (physical_page_num - 1);
    
    // Set this bit in physical bit map
	set_bit(p_bit_map, index);

    // We initialize all 2nd level page tables information in page table directtory to -1.
	struct page *ptr = &physical_memory[index];
	unsigned long *p = ptr->arr;

	for(int i = 0; i < (1 << outer_level_bits); i++)
    {
        // Initialize all 2nd level entries
        p[i] = -1;
    }
	
}

/*
Function responsible for allocating and setting your physical memory 
*/
void set_physical_mem() {

    //Allocate physical memory using mmap or malloc; this is the total size of
    //your memory you are simulating

    
    //HINT: Also calculate the number of physical and virtual pages and allocate
    //virtual and physical bitmaps and initialize them

    split_virtual_bits();
    int i = 0;

    physical_memory = (struct page *)malloc(MEMSIZE);

    // Allocate memory for physical bit map
    // Since 1 byte = 8 bits and each bit can be used to 
    // store info about one page.
    p_bit_map = (char *)malloc(physical_page_num/8);
    for(i = 0; i < sizeof(p_bit_map)/sizeof(p_bit_map[0]); i++)
    {
        p_bit_map[i] = '\0';
    }

    // Allocate memory for virtual bit map
    v_bit_map = (char *)malloc(virtual_page_num/8);
    for(i = 0; i < sizeof(v_bit_map)/sizeof(v_bit_map[0]); i++)
    {
        v_bit_map[i] = '\0';
    }

    intialize_PT();

}

static void set_bit(char *bitmap, int index)
{
    // Calculate the starting location where index is present.
    char *location = ((char *) bitmap) + (index / 8);
    
    // Find the exact bit to set
    char bit = 1 << (index % 8);

    // Set the bit at that index to 1.
    *location |= bit;
   
    return;
}

static int get_bit(char *bitmap, int index)
{
    // Calculate the starting location where index is present.
    char *location = ((char *) bitmap) + (index / 8);

    // Getting the bit at required index
    int ans = (int)(*location >> (index % 8)) & 0x1;
    
    return ans;
}

static void reset_bit(char *bitmap, int index)
{
    //Calculate the starting location where index is present.
    char *location = ((char *) bitmap) + (index / 8);
    
    // Find exact bit to reset
    char bit = 1 << (index % 8);

    // Reset the bit
    *location &= ~bit;
   
    return;
}

/*
 * Part 2: Add a virtual to physical page translation to the TLB.
 * Feel free to extend the function arguments or return type.
 */
int
add_TLB(int vpn, int ppn)
{

    /*Part 2 HINT: Add a virtual to physical page translation to the TLB */

    /*
        WE have chnaged this API signature to take VPN and PPN as input forease of implementation
    */

    int index = vpn % TLB_ENTRIES;

	tlb_store.arr[index][0] = vpn;
	tlb_store.arr[index][1] = ppn;

    return -1;
}


/*
 * Part 2: Check TLB for a valid translation.
 * Returns the physical page address.
 * Feel free to extend this function and change the return type.
 */
int
check_TLB(void *va) {

    /* Part 2: TLB lookup code here */

    unsigned long vpn = (unsigned int)va >> offset_bits;
	int index = vpn % TLB_ENTRIES;

	if (tlb_store.arr[index][0] == vpn)
    {
		unsigned long ppn = tlb_store.arr[index][1];
    	return ppn;
	}
    else
    {
		return -1;
	}

}


/*
 * Part 2: Print TLB miss rate.
 * Feel free to extend the function arguments or return type.
 */
void
print_TLB_missrate()
{
    double miss_rate = 0;	

    /*Part 2 Code here to calculate and print the TLB miss rate*/

    miss_rate = (num_tlb_misses * 1.0)/(num_tlb_lookups * 1.0);


    fprintf(stderr, "TLB miss rate %lf \n", miss_rate);
}



/*
The function takes a virtual address and page directories starting address and
performs translation to return the physical address
*/
pte_t *translate(pde_t *pgdir, void *va) {
    /* Part 1 HINT: Get the Page directory index (1st level) Then get the
    * 2nd-level-page table index using the virtual address.  Using the page
    * directory index and page table index get the physical address.
    *
    * Part 2 HINT: Check the TLB before performing the translation. If
    * translation exists, then you can return physical address from the TLB.
    */

    unsigned int virt_add = (unsigned int)va;
	
    // calculate the offset mask and offset value
    unsigned long offset_mask = (1 << offset_bits);
    offset_mask -= 1;
    unsigned long offset = virt_add & offset_mask;

    num_tlb_lookups++;

    // TODO
    // Check if the corresponding PA is present in the tlb
	int tlb_page_num = check_TLB(va);

	if (tlb_page_num != -1)
    {
		char *tlb_pa = (char *) &physical_memory[tlb_page_num];
		tlb_pa = tlb_pa + offset;
		return (pte_t *)tlb_pa;
	}

    // Update that its a TLB miss
    num_tlb_misses++;

    // Calculate the outer index
    unsigned int num_bits = (32 - outer_level_bits);
    unsigned long outer_index = (virt_add >> num_bits);

    // Calculate the vpn so we can retrieve inner index
    unsigned int vpn = (virt_add >> offset_bits);
    unsigned long inner_bits_mask = (1 << inner_level_bits);
    inner_bits_mask -= 1;

    // Computing the inner index
    unsigned long inner_index = vpn & inner_bits_mask;

    // The actual index, ie which index where it is located in physical memory
    int index = (outer_index * (1 << inner_level_bits)) + inner_index;

    // Check if this index is valid using bitmap
	int bit = get_bit(v_bit_map, index);
    if (bit != 1)
    {
        return NULL;
    }

    // Calculate the 2nd level page table entry
    pde_t *pd_entry = (pgdir + outer_index);

    // The pd_entry has the corresponding inner level page number
    unsigned long inner_level_page_num = *pd_entry;

    // Compute address of inner page table
    pte_t *inner_page_table_address = (pte_t *)&physical_memory[inner_level_page_num];
    
    pte_t *pt_entry = (inner_page_table_address + inner_index);

    // Get the actual page that resides in the physical memory
    unsigned long num = *pt_entry;

    // Get the physical page address of the above page
    pte_t *physical_page_address = (pte_t *)&physical_memory[num];

    // Finally, the physical page address + offset gives the corresponding physical address for passed virtual address
    unsigned long physical_address = (unsigned long)((char *)physical_page_address + offset);
	
	// Add to the TLB
	add_TLB(vpn, num);

    return (pte_t *)physical_address;

    //If translation not successful, then return NULL
    //return NULL; 
}


/*
The function takes a page directory address, virtual address, physical address
as an argument, and sets a page table entry. This function will walk the page
directory to see if there is an existing mapping for a virtual address. If the
virtual address is not present, then a new entry will be added
*/
int
page_map(pde_t *pgdir, void *va, void *pa)
{

    /*HINT: Similar to translate(), find the page directory (1st level)
    and page table (2nd-level) indices. If no mapping exists, set the
    virtual to physical mapping */
    
    unsigned int virt_add = (unsigned int)va;

    // calculate the offset mask and offset value
    unsigned long offset_mask = (1 << offset_bits);
    offset_mask -= 1;
    unsigned long offset = virt_add & offset_mask;

    // Calculate the ppn and vpn
    unsigned long ppn = (unsigned int)pa >> offset_bits;
    unsigned long vpn = virt_add >> offset_bits;

    // Calculate the outer index
    unsigned int num_bits = (32 - outer_level_bits);
    unsigned long outer_index = (virt_add >> num_bits);

    unsigned long inner_bits_mask = (1 << inner_level_bits);
    inner_bits_mask -= 1;
    // Computing the inner index
    unsigned long inner_index = vpn & inner_bits_mask;
    

	//check if this virtual to physical page translation is already in the TLB

	int check_page = check_TLB(va);
	num_tlb_lookups++;

    // Already present in tlb
	if(check_page == ppn)
    {
        return 0;
    }

    // Increment the tlb miss since we dint find it in TLB
	num_tlb_misses++;

    // Calculate the 2nd level page table entry
    pde_t *pd_entry = (pgdir + outer_index);

    // If no mapping ,find a page moving from last but 1 page in memory
    // this will be used for 2nd level page table.
    if (*pd_entry == -1)
    {
		int last_page = physical_page_num - 1;
		while(last_page >= 0)
        {
			// Use bit map to know if there is an empty page
			int bit = get_bit(p_bit_map, last_page);
			if(bit == 0)
            {
                // Found the page!! Mark the page.
				set_bit(p_bit_map, last_page);
				// Set this page as the 2nd level page directory
				*pd_entry = last_page;
				break;
			}
			last_page--;
		}
	}

    // The pd_entry has the corresponding inner level page number
    unsigned long inner_level_page_num = *pd_entry;

    // Compute address of inner page table
    pte_t *inner_page_table_address = (pte_t *)&physical_memory[inner_level_page_num];

    pte_t *pt_entry = (inner_page_table_address + inner_index);
		
    // Store the ppn into the page table entry.
    *pt_entry = ppn;

    // Add the vpn to ppn translation to the 
	add_TLB(vpn, ppn);	
    	return 1;

    //return -1;
}


/*Function that gets the next available page
*/

//** This function return type is changed to unsigned long
//   which returns the first free page.
unsigned long get_next_avail(int num_pages) {
 
    //Use virtual address bitmap to find the next free page

    int page_start = 0;
    int i = 0;

    // A loop that finds out if there are continous free pages
    // We find the first free page, and then look up for remaining free
    // pages, if we are able to find them, we return the first index.
    while(i < virtual_page_num)
    {
        int bit = get_bit(v_bit_map, i);
        if(bit == 0)
        {
            int count = 1;
            int j = i + 1;

            while(j < virtual_page_num && count < num_pages)
            {
                bit = get_bit(v_bit_map, j);
                if(bit == 1)
                {
                    break;
                }
                else
                {
                    count++;
                    j++;
                }
            }
            if(count == num_pages)
            {
                // This means we found continous virtual pages for required memory,
                // return the first page index!
                page_start = i;
                return page_start;
            }
            i = j;
            continue;
        }
        i++;
    }
    return -1;
}


/* Function responsible for allocating pages
and used by the benchmark
*/
void *t_malloc(unsigned int num_bytes) {

    /* 
     * HINT: If the physical memory is not yet initialized, then allocate and initialize.
     */

   /* 
    * HINT: If the page directory is not initialized, then initialize the
    * page directory. Next, using get_next_avail(), check if there are free pages. If
    * free pages are available, set the bitmaps and map a new page. Note, you will 
    * have to mark which physical pages are used. 
    */

   	pthread_mutex_lock(&mutex);
	
	if (init == false)
    {
		set_physical_mem();
		init = true;
	}
	
    // Get number of pages that correspond to num_bytes passed to t_malloc
    int pages = num_bytes / (PGSIZE);
    int rest = num_bytes % (PGSIZE);
    
    // Accomodate one more page for remaining bytes
    if( rest > 0)
    {
        pages++;
    }

    // We create an array to store all the physical pages that are free that can be allocated
    int physical_mem_pages[pages];
    int cnt = 0;
    int i = 0;

    // Find the free pages in the physical memory
    while(cnt < pages && i < physical_page_num)
    {
        int bit = get_bit(p_bit_map, i);
        if(bit == 0)
        {
            physical_mem_pages[cnt] = i;
            cnt++;
        }
        i++;
    }

    // We have failed to find free physical pages that are required to alloacte a given memory
    if(cnt < pages)
    {
        pthread_mutex_unlock(&mutex);
        return NULL;
    }

    // Find the corresponding first page in virtual memory. We should be able to 
    // find a continous "pages" in virtual space.
    int start_page = get_next_avail(pages);

    // Failed to find virtual pages
    if(start_page == -1)
    {
        pthread_mutex_unlock(&mutex);
        return NULL;
    }

    int end_page = start_page + pages - 1;

    cnt = 0;

    // Find the corresponding virtual and physical address and use page_map
    // API to map these both.
    for(i=start_page; i<= end_page; i++)
    {
        unsigned long temp_virt_addr = i << offset_bits;
        unsigned long temp_phy_addr = physical_mem_pages[cnt] << offset_bits;

        set_bit(v_bit_map, i);

        set_bit(p_bit_map, physical_mem_pages[cnt]);

        cnt++;

        page_map((pde_t *)(physical_memory+physical_page_num-1), (void *)temp_virt_addr, (void *)temp_phy_addr);
    }

    // Return the virtual address ie vpn + offset
    // The vpn is the first page we found in virtual memory
    // to allocate the required memory.
    unsigned long VA = start_page << offset_bits;

    // Release the mutex
    pthread_mutex_unlock(&mutex);

    return (void *)VA;

}

/* Responsible for releasing one or more memory pages using virtual address (va)
*/
void t_free(void *va, int size) {

    /* Part 1: Free the page table entries starting from this virtual address
     * (va). Also mark the pages free in the bitmap. Perform free only if the 
     * memory from "va" to va+size is valid.
     *
     * Part 2: Also, remove the translation from the TLB
     */

    pthread_mutex_lock(&mutex);
	
    // Get number of pages to free
    int pages = size / (PGSIZE);
    int rest = size % (PGSIZE);
    
    // Accomodate one more page for remaining bytes
    if(rest > 0)
    {
        pages++;
    }

    // Get the start page from the virtual address
    unsigned long start_page = (int)va >> offset_bits;

    bool valid = true;

    for(int i = start_page; i < (start_page + pages); i++)
    {
        // Get the bit from the virtual bit map
        int bit = get_bit(v_bit_map, i);
        if(bit == 0)
        {
            valid = false;
            break;
        }
    }

    // Return if it is not valid
    if(valid == false)
    {
        pthread_mutex_unlock(&mutex);
        return;
    }

    // compute the corresponding physical and virtual address and 
    // mark those pages to freei ie reset the corresponding bit in 
    // both physical and virtual bit maps.
    for(int i = start_page; i < (start_page + pages); i++)
    {
        void * virt_addr = (void *)(start_page << offset_bits);

        pte_t phy_addr = (pte_t)(translate((pde_t *)(physical_memory + physical_page_num - 1), virt_addr));

        unsigned long physical_page = (phy_addr >> offset_bits);

        reset_bit(v_bit_map, i);

        reset_bit(p_bit_map, physical_page);

    }

    // Remove the translation from the TLB if exists.
    for(int i = start_page; i < (start_page + pages); i++)
    {
		int index = i % TLB_ENTRIES;
		if(tlb_store.arr[index][0] == i)
        {
			tlb_store.arr[index][0] = -1;
			tlb_store.arr[index][1] = -1;
		}
		
	}

    // Release the mutex
    pthread_mutex_unlock(&mutex);
}


/* The function copies data pointed by "val" to physical
 * memory pages using virtual address (va)
 * The function returns 0 if the put is successfull and -1 otherwise.
*/
void put_value(void *va, void *val, int size) {

    /* HINT: Using the virtual address and translate(), find the physical page. Copy
     * the contents of "val" to a physical page. NOTE: The "size" value can be larger 
     * than one page. Therefore, you may have to find multiple pages using translate()
     * function.
     */

    pthread_mutex_lock(&mutex);
    char *value = (char *)val;
    //Compute the physical address from virtual address using translate API
    pte_t pa = (pte_t)(translate((pde_t *)(physical_memory + physical_page_num - 1), va));
    //We will be putting the data byte wise, hence take below variables for that purpose
    char *phy_addr = (char *)pa;
    char *virt_addr = (char *) va;	
    	
	
	char *last_virt_addr = virt_addr + size;
	//checking for the validity of the virtual pages from va to va+size

	// Compute the corresponding VPN's
	unsigned int vpn_first = (unsigned int)virt_addr >> offset_bits;
	unsigned int vpn_last = (unsigned int)last_virt_addr >> offset_bits;

    // Check if these are valid VPN's
	for(int i = vpn_first; i <= vpn_last; i++)
    {
        // Get the bit from the virtual bit map
		int bit = get_bit(v_bit_map, i);

        // This means its not valid and hence return
		if (bit == 0)
        {   
			pthread_mutex_unlock(&mutex);
			return;
		}
	}

    // Since all pages are valid, we start putting the value byte wise

    int temp = 0;
    while(temp < size)
    {
        //Copy the value byte wise. This logic is similar to memcpy implementaion!!
        *phy_addr = *value;
        //Increment the pointers for remaining copy
        virt_addr++;
        phy_addr++;
        value++;
        temp++;

        unsigned int vai = (unsigned int)virt_addr;
		
        int outer_bits_mask = (1 << offset_bits);
        outer_bits_mask -= 1;

        int offset = vai & outer_bits_mask;

        // if offset = 0 then we need to update the physical page by calling translate API again
        if (offset == 0)
        {
			//update the physical page
            pa = (pte_t)(translate((pde_t *)(physical_memory + physical_page_num - 1), (void *) virt_addr));
			phy_addr = (char *) pa;
        }
    }

	pthread_mutex_unlock(&mutex);

}


/*Given a virtual address, this function copies the contents of the page to val*/
void get_value(void *va, void *val, int size) {

    /* HINT: put the values pointed to by "va" inside the physical memory at given
    * "val" address. Assume you can access "val" directly by derefencing them.
    */

   	pthread_mutex_lock(&mutex);

    char *value = (char *)val;
    //Compute the physical address from virtual address using translate API
    char * pa = (char *)translate((pde_t *)(physical_memory + physical_page_num - 1), va);
    //We will be putting the data byte wise, hence take below variables for that purpose
    char *virt_addr = (char *) va;	
    	
	
	char *last_virt_addr = virt_addr + size;
	//checking for the validity of the virtual pages from va to va+size

	// Compute the corresponding VPN's
	unsigned int vpn_first = (unsigned int)virt_addr >> offset_bits;
	unsigned int vpn_last = (unsigned int)last_virt_addr >> offset_bits;


    // Check if these are valid VPN's
	for(int i = vpn_first; i <= vpn_last; i++)
    {
        // Get the bit from the virtual bit map
		int bit = get_bit(v_bit_map, i);

        // This means its not valid and hence return
		if (bit == 0)
        {   
			pthread_mutex_unlock(&mutex);
			return;
		}
	}
    
    for (int i = 0; i < size; i++)
    {
        //Copy the value from physical address byte wise
		*value = *pa;
        // Increment the pointer for remaining copy
		value++;
		pa++;
		virt_addr++;

		
        unsigned int vai = (unsigned int)virt_addr;
		
        int outer_bits_mask = (1 << offset_bits);
        outer_bits_mask -= 1;

        int offset = vai & outer_bits_mask;

        // if offset = 0 then we need to update the physical page by calling translate API again
        if (offset == 0)
        {
			//update the physical page
            pa = (char *)(translate((pde_t *)(physical_memory + physical_page_num - 1), (void *) virt_addr));
        }

	}

	pthread_mutex_unlock(&mutex);

}



/*
This function receives two matrices mat1 and mat2 as an argument with size
argument representing the number of rows and columns. After performing matrix
multiplication, copy the result to answer.
*/
void mat_mult(void *mat1, void *mat2, int size, void *answer) {

    /* Hint: You will index as [i * size + j] where  "i, j" are the indices of the
     * matrix accessed. Similar to the code in test.c, you will use get_value() to
     * load each element and perform multiplication. Take a look at test.c! In addition to 
     * getting the values from two matrices, you will perform multiplication and 
     * store the result to the "answer array"
     */
    int x, y, val_size = sizeof(int);
    int i, j, k;
    for (i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            unsigned int a, b, c = 0;
            for (k = 0; k < size; k++) {
                int address_a = (unsigned int)mat1 + ((i * size * sizeof(int))) + (k * sizeof(int));
                int address_b = (unsigned int)mat2 + ((k * size * sizeof(int))) + (j * sizeof(int));
                get_value( (void *)address_a, &a, sizeof(int));
                get_value( (void *)address_b, &b, sizeof(int));
                // printf("Values at the index: %d, %d, %d, %d, %d\n", 
                //     a, b, size, (i * size + k), (k * size + j));
                c += (a * b);
            }
            int address_c = (unsigned int)answer + ((i * size * sizeof(int))) + (j * sizeof(int));
            // printf("This is the c: %d, address: %x!\n", c, address_c);
            put_value((void *)address_c, (void *)&c, sizeof(int));
        }
    }
}



